# To do

Planned work, roughly in order.

- **Library use.** Make the crate a first-class library, not just the
  binary's internals. The surface is now settled and narrowed (2026-09-23,
  below); what remains is to keep it stable across releases, document every
  public item with an example, and decide whether to offer C and Python
  bindings.
- [x] 2026-10-05: Wrote a code example for every public item. 97 new doctests cover all
  91 genuinely public items; 192 tests pass in total. The earlier "roughly 70" and "258"
  counts were both wrong, because they included `pub` items inside private modules
  (`distance::dna`, `distance::gaps`, `distance::bl45`) and `pub fn` inside `pub(crate)`
  structs, none of which are reachable from outside the crate.
- [x] 2026-10-05: Settled `SubstitutionMatrix`'s three failure conventions. They stay
  as they are (`blosum62`/`blosum45` panic through `.expect`, `by_name` returns
  `Option`, `from_file`/`parse` return `Result`), which Travis agreed to; the compiled-in
  matrices cannot fail, and `by_name`'s `None` is a lookup answer rather than an error.
  The caller was the problem: see the `--matrix` seam fix below.
- The minimum Rust version (1.85) is set by clap 4.6; CI builds with it.
- [x] 2026-11-06: Sorted per-join batches, packed to eight bytes an entry with the new
  node stored once per run. 20,000 taxa, load about 2, `-T 8`, baseline and new
  alternated: join loop 16.9 s to 14.6 s, peak 4.22 GB to 3.61 GB, identical tree.
  Phases after, printed under `-v`: pull 1.95 s, update 3.99 s, insert 3.35 s (bucketing
  1.06 s), rebuild 1.87 s. The sorted batches alone, with 12-byte entries, gained nothing:
  the recorded "4 s of sift-up" was wrong, since a sift-up stops at the first parent that
  is not larger and new distances are mid-range, so pushes were already near O(1). The
  cost is writing 2e8 entries per run, which packing cuts by a third. Microbenchmark of
  the sort alone, 2e8 entries in batches of 700: 12-byte copy and stable sort 3.84 s, in
  place unstable 3.11 s, packed u64 2.34 s, no sort 0.43 s.
- [x] 2026-11-06, decided against: a square matrix for the update loop. Measured the
  loop per element with the matrix in L3 and out of it: 6,000 taxa (72 MB) 12.9 ns,
  12,000 (288 MB) 15.0 ns, 20,000 (800 MB) 18.2 ns. The memory stall is about 5 of 18 ns,
  so streaming reads could recover at most 1.0 s of 14.6, and the symmetric write would
  become a scattered store and take part of that back; 2-5% for doubling the matrix.
- [x] 2026-11-06: Hoisted the update loop onto local slices with a row-base table for
  the matrix index. The disassembly showed the real cost: the orientation swap was already
  a `cmov`, but every iteration reloaded `self`'s fields from the stack and wrote the
  `fresh` length back to memory, with an `imul` per matrix access. Branch `update-loop`,
  alternated against main at load 45: update 3.67 s to 3.10 s (18.2 to 15.5 ns per
  element), join loop 14.15 s to 13.35 s, tree byte-identical.
- [x] 2026-11-06: Fused the bucketing pass into the update loop. The `fresh` buffer was
  written at 12 bytes an entry and read back by a second pass; now each new distance is
  packed straight into a bucket for the other node's cluster, and the queue mapping is
  applied after the loop once the new node's cluster is known. Alternated at load 58:
  insert 3.33 s to 2.30 s, update 3.15 s to 3.65 s (the bucket push costs more than a
  stream push), join loop 13.5 s to 13.0 s, byte-identical tree. The earlier note blaming
  the per-batch allocation was wrong: 30 small allocations a join is 0.06 s.
- In-memory engine speed, done for now: 16.9 s to 13.0 s at 20,000 taxa on 2026-11-06.
  What is left is under a second (bucket pushes into 30 buffers, bounds checks, four
  dependent loads an element), and I think the next lever is algorithmic.
- Reduce in-memory engine memory: the queues hold an entry for every pair
  (about 12 bytes each); 104 GB at 100,000 taxa.
- [x] 2026-11-06: Measured the external-memory engine with heap timers (branch
  `extmem-merges`), 20,000 taxa at `--memory 0.1`, 83.3 s at load 60: level merges 44.0 s
  in 1,744 calls, of which 9.8 s is 2.8 M one-block reads and 268 k writes; spill sorts
  4.1 s; stores 2.7 s; engine phases pull 8.4 s, update 8.1 s, insert 31.0 s, rebuild
  7.3 s; about 28 s outside the loop in the initial fill. 15.1 GB written for 4.7e8
  entries inserted, so each entry is written about 2.7 times. Merge CPU was 27 ns an
  entry: nested `Vec<Vec<_>>` bookkeeping indexed several times an entry, two random
  loads into `redirect` for the dead check, and a pop plus a push on the source heap.
- [x] 2026-11-06: Rewrote `merge_levels` over local cursors with `replace_top`: merges
  44.0 s to 31.4 s, run 83.3 s to 71.1 s, reproducible to 0.1 s across alternated pairs,
  byte-identical tree, same I/O and memory.
- [x] 2026-11-06: Positioned reads and writes (`FileExt::read_at`/`write_all_at`) into a
  reused byte buffer for the heap's scratch file. I/O 9.8 s to 9.0 s, run 71.1 s to
  70.1 s: the seek was cheap, and what remains per 4 KB read is the read itself. The first
  version kept the buffer at its largest size, a whole run from `store`, which added
  4 MB of peak memory across 28 heaps that the plan did not account for, and bounding it
  to a merge's output chunk still left 1.5 MB; it is now kept only at block size.
- [x] 2026-11-06, decided against: a bitset of live nodes for the dead check in merges,
  instead of two loads into `redirect`. Alternated against the committed binary: merges
  31.95/32.12 s without, 32.39/32.85 s with. Out-of-order execution already hides those
  loads. Reverted; the merge loop's remaining cost is the source-heap sift and the copy.
- [x] 2026-11-06: Measured the pager/flush change (64 KB row pagers reused across joins,
  flush walking the `rows` array, `PAGER_BYTES` in the plan), two alternated pairs at
  load 40 on a noisy machine: matrix reads 3.8 M to 1.38 M (4.0 s to 2.6-3.0 s), flushes
  3.5 s to 2.7-3.4 s, update phase 12.8/9.1 s to 7.8/9.9 s, run 79.9/73.8 s to
  71.8/76.0 s, byte-identical tree, same peak. Signal about 2-3 s, noise about 2 s.
  Kept: memory-trivial and accounted. Committed as `ea6676b` on `extmem-merges`.
- [x] 2026-11-06: Priced 30% more memory for the read-chunk decision: the committed
  binary at `--memory 0.13` ran 69.3 s (peak 90.5 MB) against 73.8-79.9 s at 0.1, but
  merges were 31.9 s against 32.2-32.8 s, unchanged; the gain was in the pull phase,
  8.5 s to 5.0 s, from the plan keeping more clusters. So spending ~30% of a heap's
  allowance on merge read chunks costs well under a second of merge time while the
  chunks could save 3-5 s of the 10 s of heap I/O. Read chunks are worth building,
  with the budget taught about them. That the pull phase is that sensitive to the
  cluster count is a separate lever worth a look (`TARGET_STRUCT_BYTES` in budget.rs).
- [x] 2026-11-06: The rebuild's pager in `cluster_and_heap` paged at one window width,
  1.7 KB, and was allocated per row; it now borrows `self.pager_i` (64 KB). Two
  alternated pairs at load 42: matrix reads 1.38 M to 175 k, rebuild phase 6.41/8.38 s
  to 6.14/6.80 s, wall 76.0/79.1 s to 75.3/78.6 s, byte-identical tree, peak unchanged.
  Kept. The same commit splits `-v`'s heap I/O line into read and write time.
- [x] 2026-11-06: Split heap I/O time, one run at load 42: reads 3.3 s for 2.8 M
  one-block calls (1.2 us each, 11.5 GB), writes 7.8 s for 268 k ten-block calls (29 us
  each, 15.1 GB). Merge read chunks can therefore save at most about 2.5 s, below the
  noise, while a chunk of four costs about 23 MB at 20k taxa (69 live sources x 3 extra
  blocks x 28 heaps) or, paid for inside the plan, a 2x smaller run size and one more
  level. The writes are the lever: 2.7 rewrites an entry, set by the slot count.
- [x] 2026-11-06: Measured four-block merge reads (prototype, no budget change), two
  alternated pairs at load 42: read calls 2.8 M to 0.92 M but read time only 3.5 s to
  2.8-2.9 s, merges 32.6/32.3 s to 31.7/32.4 s, wall 76.1 s to 74.9/74.8 s, peak
  67.5-70.5 MB to 75.1-78.2 MB, trees identical. The read cost is the 11.5 GB copied,
  not the call count, so chunks buy about 1 s for 5-8 MB. Reverted; not worth a budget
  revamp. Dead end, with the 0.13 pricing above as its false lead.
- External-memory engine, what is left at 20,000 taxa and 70 s: merges 31 s, of which
  I/O 9 s in 2.8 M one-block reads; each entry is written about 2.7 times. More slots
  per level would cut promotions (with 17 slots at this allowance, two per entry; with
  64, one) at the cost of a head block per slot resident in `h2`, which the plan must
  account for. Engine phases pull 8.3 s, update 8.1 s (the disk-window update loop, 40 ns
  an element, not yet hoisted like the in-memory one), initial fill about 28 s.
- External-memory engine, 20,000 taxa at a 100 MB budget (85.7 s on an idle machine,
  measured 2026-10-05): level
  merges of the disk heaps dominate, then pulls, the update loop, spill
  sorts and staging pushes. More slots per level would reduce how many
  merges each entry passes through.
- [x] 2026-10-05: Clamped negative branch lengths to zero in both engines, with the
  count reported on stderr. The framing in this item was wrong: both children can never
  be negative, since their lengths sum to twice the pair distance. What happens is that
  the existing fold moves a negative length onto its sibling, and when the pair distance
  is itself negative the sibling stays negative. Reference outputs keep their negative
  values and the test helpers clamp both sides, so they still work as an independent
  oracle; see `docs/testing.md`.
- Decide whether the 100 MB minimum budget should scale with the taxon count. At 100,000
  taxa the plan reports 121 MB as the true minimum and the run sits at 132 MB, so a flat
  floor understates what large inputs need. Raised 2026-09-22.
- The budget split (7/20 window, 1/10 build scratch, 3/20 candidate structures, the rest to
  cluster-pair heaps) is a guess, not a measured optimum. Raised 2026-09-22.
- [x] 2026-09-24: Measured what a tight budget costs. Against a roomy configuration it is
  about 1.4x, flat with size: 1.45x at 20,000 taxa and 1.38x at 50,000. The cluster count
  is most of it (4 to 30 clusters buys 26%); window width buys 13% at 20,000 taxa, and
  4096 columns was no better than 1024 at 50,000. Keeping all 30 clusters at 100,000 taxa
  would need about 1.24 GB against today's 121 MB, so roughly 10x the memory for 1.4x the
  speed. Travis decided to keep the floor at feasibility. Note the caveat: this machine
  holds the whole matrix file in page cache, so a window miss costs a syscall rather than a
  seek; on a memory-tight machine the window would matter more, and dropping the page cache
  needs root.
- Consider a verbose line saying how much more memory would buy, so the 1.4x is visible
  rather than hidden. Raised 2026-09-24.
- [x] 2026-09-25: Cut the plan's under-accounting from 5.8% to 1.7%. At 100,000 taxa it
  now raises to 124 MB and the run peaks at 126.1 MB (was 121 planned, 128 actual).
  Added the per-heap staging buffers, the join loop's per-taxon scratch (`rows`,
  `horiz`), and the true candidate slot cost (13 to 29 bytes, covering `free_cands` and
  the drain buffer); `ArrayHeap::run_size_for` sizes staging exactly, and the cluster
  reduction now targets `min_struct + min_stage`.
- Decide whether the last 1.7% is worth chasing. I think it is allocator slack and the
  leaf-name strings, which the plan neither allocates nor controls, and recommended
  closing it; Travis has not ruled. Raised 2026-09-25.
- [x] 2026-09-25: Fixed the `--matrix` seam. An unknown built-in name fell through to
  the file reader, so `--matrix BLOSUM80` complained about a missing file. It now
  reports the unknown name and lists the built-ins; `SubstitutionMatrix::names()` keeps
  the message and the `match` in step. The three failure conventions themselves stay as
  they are, which Travis agreed to.
- [x] 2026-10-06: Released 2.1.0, tagged `v2.1.0`, on crates.io, docs.rs built, GitHub
  release published. 2.0.0 stays unyanked: it had 0 downloads, but crates.io cannot
  replace a published version and resolution picks 2.1.0 regardless.
- [x] 2026-10-05: Released 2.0.0, tagged `v2.0.0` and on crates.io. Made `nj::extmem`
  internal first, so `MemoryPlan`'s fields are not frozen by the stable API and the
  accounting can still change in 2.x.
- [x] 2026-09-25: Released 2.0.0-rc.4, tagged and on crates.io.
- [x] 2026-09-24: A 100,000-taxon run completes: 51m56s, 128 MB peak, 99,999 joins, a tree
  with 100,000 distinct labels. It had never finished before; see the freeze hang below.
- [x] 2026-09-24: Fixed the hang in the candidate-heap drain (appending from inside the
  loop over `cand_heaps` re-indexed it and, with one heap allowed, fed pairs straight back
  into the heap they came from). Regression test added; verified it hangs without the fix.
- Benchmark and probe runs must cap worker threads (`ninja -T N`, or a rayon pool in a test
  binary). The distance step otherwise takes all 192 cores on this shared machine, which
  disrupts other users. Raised 2026-09-24.
- [x] 2026-10-05: Got the 20,000-taxon time on an idle machine (load 0.4): 85.7 s,
  70.1 MB peak, `--memory 0.1 -T 8`. Earlier figures ranged 80 s to 445 s purely on
  other users' load, so this is the one to quote.
- Find a way to run multi-hour jobs that the harness watchdog will not kill. It killed the
  100,000-taxon run twice for "low memory" while the machine had 572 GB available. A
  `setsid nohup` relaunch on 2026-09-23 is the current attempt; tmux and systemd-run are
  also on the machine and untried. Raised 2026-09-23.
- [x] 2026-09-23: Narrowed the public API surface, from 258 items to about
  70. `heap`, `distance::gaps`, `nj::extmem::matrix` and `nj::extmem::budget`
  are private modules; the `CandidateHeap` re-export is gone; `cluster::write_table`,
  `DistanceMatrix::prefetch` and the raw `DiskMatrix` accessors and layout
  constants are `pub(crate)`. `NjStats` is re-exported at the crate root.
  `DiskMatrix::prefetch`, `DiskMatrix::mem_row`, `MinHeap::iter` and
  `CandidateHeap::is_empty` turned out to be uncalled and were deleted;
  four accessors used only by unit tests are now `#[cfg(test)]`. Not yet committed.
- [x] 2026-09-23: Deleting the uncalled `DiskMatrix::prefetch` removed the
  crate's second `unsafe` block, so the claim at `src/lib.rs:45` that there
  is only one is true again.
- [x] 2026-10-05: Released 2.0.0 final. Travis first asked on 2026-09-21 and pulled it
  back when he saw the open items; the API review finished with the rc.4 narrowing and
  making `nj::extmem` internal.

- [x] 2026-09-15: `-v` with no value fails ("a value is required for --verbose"); accept a bare `-v` as verbose level 2 (clap default_missing_value). Seen by a colleague 2026-09-14.
- [x] 2026-09-16, decided against: Quote Newick labels that contain `#` (and any other character outside the plain set): extended-Newick readers such as IcyTree and Dendroscope read `name#tag` as a reticulation node and merge every leaf sharing the tag, so RepeatMasker-style names come out as a network with cycles. Found 2026-09-14. Quoting only helps IcyTree; Dendroscope and SplitsTree strip quotes and still merge labels whose tag starts with H, L or R, so Travis dropped the idea as a hack with little benefit.
- [x] 2026-09-21: `--collapse_identical` is now the default, with `--no_collapse_identical` to turn it off.
- [x] 2026-09-22: `--memory` now bounds the external-memory engine (20,000 taxa: 309 MB
  down to 76 MB).
