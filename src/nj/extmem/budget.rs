//! Dividing the memory budget among the external-memory engine's structures.
//!
//! The budget used to scale three things without ever multiplying by how
//! many of them there were, so the engine ran several times over it. Here
//! the per-taxon structures that cannot be made smaller are subtracted
//! first, and what is left is divided by the actual number of structures.

use crate::heap::ArrayHeap;
use crate::nj::extmem::matrix::MIN_WINDOW_COLS;

/// Smallest budget the external-memory engine accepts. Below this there is
/// not enough left, after the per-taxon structures, to give the heaps
/// buffers worth having.
pub const MIN_MEMORY_BYTES: u64 = 100 << 20;

/// Resident bytes per taxon outside the tunable structures: the tree arena
/// and its names, the redirect table, the active list, the row sums and the
/// cluster assignments. Rounded up from the measured total.
const FIXED_PER_TAXON: u64 = 176;

/// Per-structure allowance we try to reach before spending memory on more
/// cluster pairs. Smaller heaps spill to disk more often.
const TARGET_STRUCT_BYTES: u64 = 1 << 20;

/// Fewest clusters worth keeping: below this the search bounds are too
/// loose to save any work.
const MIN_CLUSTERS: usize = 4;

/// Bytes per candidate-list entry (`cand_d`, `cand_i`, `cand_j`,
/// `cand_active`).
const CAND_ENTRY_BYTES: u64 = 13;

/// Bytes each live candidate heap needs for its per-taxon vectors,
/// independent of its buffers: `r_primes` plus the row lists.
const CAND_HEAP_PER_TAXON: u64 = 32;

/// How the budget is divided. Every field is a resident allowance in bytes
/// unless its name says otherwise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MemoryPlan {
    /// The budget actually used, after raising a request below the minimum.
    pub budget: u64,
    /// Resident window of the distance matrix.
    pub window_bytes: u64,
    /// Scratch for filling the matrix: the rows being computed in parallel
    /// whose columns are bound for disk.
    pub build_scratch: u64,
    /// Each cluster-pair heap, including its staging buffer.
    pub per_pair_heap: u64,
    /// Each candidate heap's buffers, excluding its per-taxon vectors.
    pub per_cand_heap: u64,
    /// Most candidate heaps that may be live at once.
    pub max_cand_heaps: usize,
    /// Longest candidate list before a heap takes over.
    pub cand_list_cap: usize,
    /// Cluster count, reduced from the requested one when the budget cannot
    /// give `cluster_count * (cluster_count + 1) / 2` heaps a useful size.
    pub cluster_count: usize,
    /// Whether the cluster count was reduced.
    pub clusters_reduced: bool,
    /// Whether the requested budget was raised to the minimum.
    pub budget_raised: bool,
    /// Least this many taxa can be built in, whatever the budget: the
    /// per-taxon structures, the narrowest window, the fewest clusters and
    /// one candidate heap.
    pub min_required: u64,
    /// Whether `min_required` exceeds the budget, so the run will overrun it.
    pub budget_short: bool,
}

/// Heaps needed for `cc` clusters: one per unordered pair, diagonal included.
pub fn pair_heap_count(cc: usize) -> u64 {
    (cc as u64 * (cc as u64 + 1)) / 2
}

impl MemoryPlan {
    /// Divide `memory_bytes` for `k` taxa and at most `cluster_count`
    /// clusters. `input_bytes` is what the input itself occupies while the
    /// matrix is built, which the engine cannot free and so does not get to
    /// spend. `cand_heap_cap` and `cand_list_cap` are the engine's own hard
    /// limits, which the plan may lower but never raise.
    pub fn new(
        k: usize,
        cluster_count: usize,
        memory_bytes: u64,
        input_bytes: u64,
        cand_heap_cap: usize,
        cand_list_cap: usize,
    ) -> Self {
        let budget = memory_bytes.max(MIN_MEMORY_BYTES);
        let k64 = k as u64;
        let fixed = FIXED_PER_TAXON * k64 + input_bytes;
        // Never let the per-taxon structures swallow the whole budget: the
        // engine needs working room even for a very large k under a small
        // budget, and going over is better than not running.
        let free = budget.saturating_sub(fixed).max(budget / 4);

        // One candidate heap is not optional: without it the candidate list
        // has nowhere to go when it grows too long. Reserve it, and the
        // per-taxon vectors it carries, before dividing anything.
        let min_struct = ArrayHeap::min_allowance();
        let per_cand_taxon = CAND_HEAP_PER_TAXON * k64;
        let free = free.saturating_sub(per_cand_taxon + min_struct).max(budget / 8);

        let window_bytes = free * 7 / 20;
        let build_scratch = free / 10;
        let cand_share = free * 3 / 20;
        let pair_share = free - window_bytes - build_scratch - cand_share;

        let mut cc = cluster_count.max(1);
        while cc > MIN_CLUSTERS && pair_share / pair_heap_count(cc) < TARGET_STRUCT_BYTES {
            cc -= 1;
        }
        let per_pair_heap = (pair_share / pair_heap_count(cc)).max(min_struct);

        // The reserved heap, plus as many more as the candidate share pays
        // for at the same size.
        let list_share = cand_share / 2;
        let heap_share = cand_share - list_share;
        let max_cand_heaps =
            (1 + heap_share / (per_cand_taxon + min_struct)).clamp(1, cand_heap_cap as u64) as usize;
        let per_cand_heap = min_struct.max(
            (min_struct + heap_share).saturating_sub(per_cand_taxon * (max_cand_heaps as u64 - 1))
                / max_cand_heaps as u64,
        );

        let min_required = fixed
            + 4 * MIN_WINDOW_COLS as u64 * k64
            + pair_heap_count(cc) * min_struct
            + per_cand_taxon
            + min_struct
            + 1024 * CAND_ENTRY_BYTES;

        MemoryPlan {
            budget,
            window_bytes,
            build_scratch,
            per_pair_heap,
            per_cand_heap,
            max_cand_heaps,
            cand_list_cap: cand_list_cap.min((list_share / CAND_ENTRY_BYTES) as usize).max(1024),
            cluster_count: cc,
            clusters_reduced: cc < cluster_count,
            budget_raised: memory_bytes < MIN_MEMORY_BYTES,
            min_required,
            budget_short: min_required > budget,
        }
    }

    /// Lines describing any adjustment, for `--verbose`.
    pub fn notes(&self, requested: u64, requested_clusters: usize) -> Vec<String> {
        let mut out = Vec::new();
        if self.budget_raised {
            out.push(format!(
                "warning: a memory budget of {} MB is below the {} MB minimum; using {} MB",
                requested / (1 << 20),
                MIN_MEMORY_BYTES / (1 << 20),
                self.budget / (1 << 20)
            ));
        }
        if self.budget_short {
            out.push(format!(
                "warning: this many taxa need at least {} MB for the tree, the index structures \
                 and the narrowest usable matrix window; the run will exceed the {} MB budget",
                self.min_required.div_ceil(1 << 20),
                self.budget / (1 << 20)
            ));
        }
        if self.clusters_reduced {
            out.push(format!(
                "Reduced the cluster count from {} to {} so that each of the {} cluster-pair \
                 heaps fits the budget",
                requested_clusters,
                self.cluster_count,
                pair_heap_count(self.cluster_count)
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plan(k: usize, gb: f64) -> MemoryPlan {
        plan_with_input(k, gb, 0)
    }

    fn plan_with_input(k: usize, gb: f64, input: u64) -> MemoryPlan {
        MemoryPlan::new(k, 30, (gb * (1u64 << 30) as f64) as u64, input, 100, 2_000_000)
    }

    /// The parts we hand out must add up to no more than the budget.
    fn accounted(p: &MemoryPlan, k: usize) -> u64 {
        accounted_with_input(p, k, 0)
    }

    fn accounted_with_input(p: &MemoryPlan, k: usize, input: u64) -> u64 {
        input
            + FIXED_PER_TAXON * k as u64
            + p.window_bytes
            + p.build_scratch
            + p.per_pair_heap * pair_heap_count(p.cluster_count)
            + (p.per_cand_heap + CAND_HEAP_PER_TAXON * k as u64) * p.max_cand_heaps as u64
            + p.cand_list_cap as u64 * CAND_ENTRY_BYTES
    }

    #[test]
    fn allowances_stay_inside_the_budget() {
        for k in [1_000, 20_000, 100_000, 1_000_000] {
            for gb in [0.1, 0.5, 2.0, 10.0] {
                let p = plan(k, gb);
                if p.budget_short {
                    // Too many taxa for the budget whatever we do; the plan
                    // says so instead of pretending to fit.
                    assert!(p.min_required > p.budget, "k={k} gb={gb}");
                    assert!(p.notes(0, 30).iter().any(|n| n.contains("will exceed")), "k={k}");
                    continue;
                }
                assert!(accounted(&p, k) <= p.budget, "k={k} gb={gb}: {} > {}", accounted(&p, k), p.budget);
            }
        }
    }

    /// What the input occupies while the matrix is built is spent before
    /// anything else.
    #[test]
    fn the_input_comes_off_the_top() {
        let k = 20_000;
        let input = 40 << 20;
        let p = plan_with_input(k, 0.2, input);
        assert!(!p.budget_short);
        assert!(accounted_with_input(&p, k, input) <= p.budget);
        assert!(p.window_bytes < plan(k, 0.2).window_bytes, "a big input leaves a smaller window");
    }

    #[test]
    fn a_million_taxa_do_not_fit_in_the_minimum_budget() {
        let p = plan(1_000_000, 0.1);
        assert!(p.budget_short);
        // The tree and index structures alone are over 170 MB.
        assert!(p.min_required > 170 << 20, "{}", p.min_required);
        assert!(!plan(1_000_000, 10.0).budget_short);
    }

    #[test]
    fn a_small_request_is_raised_to_the_minimum() {
        let p = plan(20_000, 0.01);
        assert!(p.budget_raised);
        assert_eq!(p.budget, MIN_MEMORY_BYTES);
        assert!(!plan(20_000, 2.0).budget_raised);
        let notes = p.notes((0.01 * (1u64 << 30) as f64) as u64, 30);
        assert!(notes[0].contains("below the 100 MB minimum"), "{notes:?}");
    }

    #[test]
    fn a_tight_budget_buys_fewer_clusters() {
        let tight = plan(20_000, 0.1);
        let roomy = plan(20_000, 10.0);
        assert!(tight.cluster_count < roomy.cluster_count, "{tight:?} {roomy:?}");
        assert!(tight.cluster_count >= MIN_CLUSTERS);
        assert_eq!(roomy.cluster_count, 30, "a roomy budget keeps every cluster asked for");
        assert!(tight.clusters_reduced && !roomy.clusters_reduced);
        assert!(tight.notes(0, 30)[0].contains("Reduced the cluster count from 30"));
    }

    #[test]
    fn every_allowance_stays_usable() {
        for k in [2, 1_000, 5_000_000] {
            let p = plan(k, 0.1);
            assert!(p.per_pair_heap >= ArrayHeap::min_allowance());
            assert!(p.per_cand_heap >= ArrayHeap::min_allowance());
            assert!(p.max_cand_heaps >= 1);
            assert!(p.cand_list_cap >= 1024);
            assert!(p.window_bytes > 0);
            assert!(p.build_scratch > 0);
        }
    }
}
