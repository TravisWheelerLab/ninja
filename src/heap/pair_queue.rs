//! Priority queue of node pairs used by the in-memory engine.
//!
//! Entries are consumed in increasing distance order with lazy deletion:
//! a pair whose node has been joined is discarded when it reaches the
//! front. Most entries never move after they are inserted, so bulk inserts
//! are kept as sorted runs read through cursors: the run a rebuild
//! supplies, and one run per join holding the new node's distances. A
//! small heap keyed by each run's front merges the runs, so a pop costs a
//! cursor increment and a sift in a heap of runs rather than a sift in a
//! heap of entries, and an insert costs its share of one sort rather than
//! a sift-up. Only single inserts go into a binary heap of entries.
//!
//! A join's run stores the new node once and packs each entry into a
//! `u64`, which makes the batch two thirds the size of full pairs and lets
//! it sort as plain integers. The rebuild's run keeps full pairs, since
//! both nodes vary there.
//!
//! Ties resolve oldest first: an older run before a newer one, and any run
//! before the entry heap; within a join's run, by the other node's index.
//! In reference-order mode every entry goes through the entry heap, which
//! reproduces the original implementation's order among equal keys.

use super::MinHeap;

/// Heap key: fixed-point distance; payload: node indices `(i, j)`.
pub type Entry = (i32, (i32, i32));

/// Flipping the sign bit maps `i32` order onto `u32` order.
const SIGN: u32 = 0x8000_0000;

/// Drained batch buffers kept per queue for reuse.
const SPARE_MAX: usize = 4;

/// A sorted run: the pairs a rebuild supplies, or one join's entries
/// packed against the node they share.
#[derive(Debug, Clone)]
enum Run {
    /// Full pairs, sorted by distance.
    Pairs(Vec<Entry>),
    /// One join's entries. `node` is common to all of them; each word
    /// packs the distance and the other node, see [`PairQueue::pack`].
    Packed { node: i32, keys: Vec<u64> },
}

impl Run {
    #[inline]
    fn len(&self) -> usize {
        match self {
            Run::Pairs(v) => v.len(),
            Run::Packed { keys, .. } => keys.len(),
        }
    }

    #[inline]
    fn key(&self, i: usize) -> i32 {
        match self {
            Run::Pairs(v) => v[i].0,
            Run::Packed { keys, .. } => PairQueue::unpack_key(keys[i]),
        }
    }

    #[inline]
    fn get(&self, i: usize) -> Entry {
        match self {
            Run::Pairs(v) => v[i],
            Run::Packed { node, keys } => PairQueue::unpack(keys[i], *node),
        }
    }
}

/// A queue of pairs ordered by distance.
#[derive(Debug, Clone)]
pub struct PairQueue {
    /// Sorted runs, oldest first. An exhausted run's buffer is freed or
    /// recycled, but its slot stays so run ids remain stable.
    runs: Vec<Run>,
    /// Next unread position in each run.
    cursors: Vec<usize>,
    /// The front of every unexhausted run, keyed by `(distance, run id)` so
    /// that equal distances come out oldest run first; the payload is the
    /// run id again.
    heads: MinHeap<(i32, u32), u32>,
    /// Entries inserted one at a time.
    heap: MinHeap<i32, (i32, i32)>,
    /// Buffers of drained packed runs, handed back by `push_batch` so the
    /// caller's bucketing reuses their capacity.
    spare: Vec<Vec<u64>>,
    reference_order: bool,
}

impl PairQueue {
    /// An empty queue. With `reference_order`, bulk fills also go through
    /// the heap so that ties resolve as in the reference implementation.
    pub fn new(reference_order: bool) -> Self {
        PairQueue {
            runs: Vec::new(),
            cursors: Vec::new(),
            heads: MinHeap::new(),
            heap: MinHeap::new(),
            spare: Vec::new(),
            reference_order,
        }
    }

    /// Pack one entry of a join's batch: the distance in the high half,
    /// shifted so that `u64` order is `i32` order, and the other node with
    /// its orientation bit in the low half. `other_first` records that the
    /// pair reads `(other, node)` rather than `(node, other)`.
    #[inline]
    pub fn pack(dist: i32, other: i32, other_first: bool) -> u64 {
        debug_assert!(other >= 0);
        let hi = (dist as u32 ^ SIGN) as u64;
        let lo = (((other as u32) << 1) | other_first as u32) as u64;
        (hi << 32) | lo
    }

    #[inline]
    fn unpack_key(word: u64) -> i32 {
        ((word >> 32) as u32 ^ SIGN) as i32
    }

    #[inline]
    fn unpack(word: u64, node: i32) -> Entry {
        let lo = word as u32;
        let other = (lo >> 1) as i32;
        let pair = if lo & 1 == 1 { (other, node) } else { (node, other) };
        (Self::unpack_key(word), pair)
    }

    /// Remove every entry. The runs' memory is released, apart from a few
    /// batch buffers kept for reuse; the heaps keep their allocations.
    pub fn clear(&mut self) {
        for run in std::mem::take(&mut self.runs) {
            if let Run::Packed { keys, .. } = run {
                self.recycle(keys);
            }
        }
        self.cursors.clear();
        self.heads.clear();
        self.heap.clear();
    }

    fn recycle(&mut self, mut keys: Vec<u64>) {
        if self.spare.len() < SPARE_MAX && keys.capacity() > 0 {
            keys.clear();
            self.spare.push(keys);
        }
    }

    /// Number of entries.
    #[cfg(test)]
    pub fn len(&self) -> usize {
        let in_runs: usize = self.runs.iter().zip(&self.cursors).map(|(run, &cur)| run.len() - cur).sum();
        in_runs + self.heap.len()
    }

    /// True when no entries remain.
    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Insert one entry.
    #[inline]
    pub fn push(&mut self, key: i32, pair: (i32, i32)) {
        self.heap.push(key, pair);
    }

    /// Insert many entries at once (a rebuild). In reference-order mode
    /// they are pushed one by one, in the given order; otherwise they are
    /// sorted into a run.
    pub fn fill(&mut self, mut entries: Vec<Entry>) {
        if self.reference_order {
            for (k, p) in entries {
                self.heap.push(k, p);
            }
        } else {
            entries.sort_unstable_by_key(|e| e.0);
            self.add_run(Run::Pairs(entries));
        }
    }

    /// Insert the entries one join produces, packed with [`pack`](Self::pack)
    /// against `node`, as a single sorted run. In reference-order mode they
    /// are pushed one by one, in the given order. The batch buffer is taken
    /// and replaced with an empty one, recycled from a drained run when
    /// there is one, so the caller can refill it without reallocating.
    pub fn push_batch(&mut self, node: i32, batch: &mut Vec<u64>) {
        if self.reference_order {
            for &word in batch.iter() {
                let (k, p) = Self::unpack(word, node);
                self.heap.push(k, p);
            }
            batch.clear();
        } else {
            // Most runs outlive the next rebuild rather than draining, so a
            // spare is often missing; size the replacement to the batch just
            // taken so the caller's next bucketing does not grow from zero.
            let replacement = self.spare.pop().unwrap_or_else(|| Vec::with_capacity(batch.len()));
            let mut keys = std::mem::replace(batch, replacement);
            keys.sort_unstable();
            self.add_run(Run::Packed { node, keys });
        }
    }

    /// Append an already sorted run and register its front.
    fn add_run(&mut self, run: Run) {
        if run.len() == 0 {
            if let Run::Packed { keys, .. } = run {
                self.recycle(keys);
            }
            return;
        }
        let id = self.runs.len() as u32;
        self.heads.push((run.key(0), id), id);
        self.runs.push(run);
        self.cursors.push(0);
    }

    /// The front entry of the oldest run holding the smallest distance.
    #[inline]
    fn run_front(&self) -> Option<(Entry, u32)> {
        let &((_, id), _) = self.heads.peek()?;
        let r = id as usize;
        Some((self.runs[r].get(self.cursors[r]), id))
    }

    /// Advance run `id` past its front, re-registering its new front or
    /// releasing the run once it is exhausted.
    fn advance(&mut self, id: u32) {
        let r = id as usize;
        self.cursors[r] += 1;
        let cur = self.cursors[r];
        if cur < self.runs[r].len() {
            self.heads.replace_top((self.runs[r].key(cur), id), id);
        } else {
            self.heads.pop();
            // Release the buffer and reset the cursor together, so an
            // exhausted run reads as empty from either side.
            let done = std::mem::replace(&mut self.runs[r], Run::Pairs(Vec::new()));
            if let Run::Packed { keys, .. } = done {
                self.recycle(keys);
            }
            self.cursors[r] = 0;
        }
    }

    /// The smallest entry. On a tie between a run and the heap, the run
    /// entry (the older one) comes first.
    #[inline]
    pub fn peek(&self) -> Option<Entry> {
        let from_run = self.run_front().map(|(e, _)| e);
        match (from_run, self.heap.peek()) {
            (None, None) => None,
            (Some(r), None) => Some(r),
            (None, Some(&(k, p))) => Some((k, p)),
            (Some(r), Some(&(k, p))) => {
                if r.0 <= k {
                    Some(r)
                } else {
                    Some((k, p))
                }
            }
        }
    }

    /// Remove the smallest entry.
    #[inline]
    pub fn pop(&mut self) -> Option<Entry> {
        match (self.run_front(), self.heap.peek()) {
            (None, None) => None,
            (Some((r, id)), None) => {
                self.advance(id);
                Some(r)
            }
            (None, Some(_)) => self.heap.pop(),
            (Some((r, id)), Some(&(k, _))) => {
                if r.0 <= k {
                    self.advance(id);
                    Some(r)
                } else {
                    self.heap.pop()
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_run_and_heap_in_order() {
        for reference in [false, true] {
            let mut q = PairQueue::new(reference);
            let bulk: Vec<Entry> = (0..1000).map(|i| ((i * 7919) % 613, (i, 0))).collect();
            q.fill(bulk.clone());
            for i in 0..300 {
                q.push((i * 31) % 613, (i, 1));
            }
            let mut keys: Vec<i32> =
                bulk.iter().map(|e| e.0).chain((0..300).map(|i| (i * 31) % 613)).collect();
            keys.sort();
            assert_eq!(q.len(), keys.len());
            for k in keys {
                assert_eq!(q.peek().unwrap().0, k);
                assert_eq!(q.pop().unwrap().0, k);
            }
            assert!(q.is_empty());
            assert!(q.pop().is_none());
        }
    }

    /// Packing round-trips every field, including negative distances and
    /// both orientations, and `u64` order equals `(distance, other)` order.
    #[test]
    fn pack_round_trips_and_orders() {
        let cases = [(-5, 7, true), (-5, 3, false), (0, 0, false), (12, 1 << 30, true), (i32::MAX, 9, false)];
        for &(d, o, first) in &cases {
            let (k, pair) = PairQueue::unpack(PairQueue::pack(d, o, first), 42);
            assert_eq!(k, d);
            assert_eq!(pair, if first { (o, 42) } else { (42, o) });
        }
        let mut words: Vec<u64> = cases.iter().map(|&(d, o, f)| PairQueue::pack(d, o, f)).collect();
        words.sort_unstable();
        let order: Vec<(i32, i32)> =
            words.iter().map(|&w| (PairQueue::unpack_key(w), (w as u32 >> 1) as i32)).collect();
        let mut want: Vec<(i32, i32)> = cases.iter().map(|&(d, o, _)| (d, o)).collect();
        want.sort();
        assert_eq!(order, want);
    }

    /// Many batches of varying size, interleaved with singles and pops, in
    /// both modes: every key comes out in order and nothing is lost.
    #[test]
    fn merges_many_batches_and_singles() {
        for reference in [false, true] {
            let mut q = PairQueue::new(reference);
            let mut expected: Vec<i32> = Vec::new();
            let mut popped: Vec<i32> = Vec::new();
            let mut seed: u32 = 12345;
            let mut next = move || {
                seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
                (seed >> 16) & 0x7fff
            };
            q.fill((0..500).map(|i| ((i * 7919) % 613, (i, 0))).collect());
            expected.extend((0..500).map(|i| (i * 7919) % 613));
            let mut batch: Vec<u64> = Vec::new();
            for round in 0..200 {
                let size = (next() % 40) as usize; // 0..39, so some batches are empty
                batch.clear();
                for j in 0..size {
                    let d = (next() % 613) as i32 - 300; // negatives included
                    expected.push(d);
                    batch.push(PairQueue::pack(d, j as i32, j % 2 == 0));
                }
                q.push_batch(1000 + round, &mut batch);
                assert!(batch.is_empty(), "push_batch leaves the caller an empty buffer");
                if round % 3 == 0 {
                    let k = (next() % 613) as i32 - 300;
                    expected.push(k);
                    q.push(k, (round, -1));
                }
                let mut prev = i32::MIN;
                for _ in 0..(next() % 5) {
                    if let Some((k, _)) = q.pop() {
                        assert!(k >= prev, "pops within a burst are ordered");
                        prev = k;
                        popped.push(k);
                    }
                }
            }
            let mut prev = i32::MIN;
            while let Some((k, _)) = q.pop() {
                assert!(k >= prev);
                prev = k;
                popped.push(k);
            }
            assert!(q.is_empty());
            expected.sort();
            popped.sort();
            assert_eq!(popped, expected, "a key was lost or invented");
        }
    }

    /// Equal distances come out oldest first: the rebuild run, then each
    /// batch in insertion order, then the heap; within a batch, by the
    /// other node's index.
    #[test]
    fn ties_come_out_oldest_first() {
        let mut q = PairQueue::new(false);
        q.fill(vec![(5, (0, 0)), (5, (0, 1))]);
        let mut b1 = vec![PairQueue::pack(5, 9, false), PairQueue::pack(5, 2, true)];
        q.push_batch(100, &mut b1);
        q.push(5, (3, 0));
        let mut b2 = vec![PairQueue::pack(5, 0, false)];
        q.push_batch(200, &mut b2);
        let order: Vec<(i32, i32)> = std::iter::from_fn(|| q.pop()).map(|e| e.1).collect();
        assert_eq!(order[0].0, 0, "rebuild run first");
        assert_eq!(order[1].0, 0);
        assert_eq!(order[2..4], [(2, 100), (100, 9)], "first batch, by other node, orientation kept");
        assert_eq!(order[4], (200, 0), "second batch before the heap");
        assert_eq!(order[5], (3, 0), "heap entry last");
    }

    /// A drained packed run gives its buffer back for the next batch; a
    /// drained pair run is freed; popping past the end is clean.
    #[test]
    fn drained_runs_are_recycled_or_released() {
        let mut q = PairQueue::new(false);
        let mut batch: Vec<u64> = Vec::with_capacity(64);
        batch.push(PairQueue::pack(1, 0, false));
        q.push_batch(7, &mut batch);
        assert!(batch.is_empty() && batch.capacity() >= 1, "no spare yet: a fresh buffer sized to the batch");
        q.fill(vec![(2, (1, 0)), (3, (1, 1))]);
        assert_eq!(q.pop(), Some((1, (7, 0))));
        assert_eq!(q.spare.len(), 1, "the drained packed run's buffer is kept");
        batch.push(PairQueue::pack(4, 5, true));
        q.push_batch(8, &mut batch);
        assert_eq!(batch.capacity(), 64, "the caller gets the recycled buffer back");
        assert_eq!(q.pop(), Some((2, (1, 0))));
        assert_eq!(q.pop(), Some((3, (1, 1))));
        if let Run::Pairs(v) = &q.runs[1] {
            assert!(v.is_empty() && v.capacity() == 0, "exhausted pair run freed");
        } else {
            panic!("run 1 should be the pair run");
        }
        assert_eq!(q.pop(), Some((4, (5, 8))));
        assert_eq!(q.pop(), None);
        assert_eq!(q.peek(), None);
        assert!(q.is_empty());
    }
}
