# To do

Planned work, roughly in order.

- **Library use.** Make the crate a first-class library, not just the
  binary's internals. The surface is now settled and narrowed (2026-09-23,
  below); what remains is to keep it stable across releases, document every
  public item with an example, and decide whether to offer C and Python
  bindings.
- Write a code example for each public item. All 258 items had a doc
  comment, enforced by `missing_docs`, but not one had a ```-fenced example;
  the only runnable ones are the two at the crate level in `src/lib.rs`.
  Roughly 70 items remain public after the narrowing. Raised 2026-09-23.
- Settle `SubstitutionMatrix`'s three failure conventions: `blosum62` and
  `blosum45` panic through `.expect`, `by_name` returns `Option`, and
  `from_file` and `parse` return `Result` (`src/distance/submat.rs:36-73`).
  Defensible as it stands, since the compiled-in matrices cannot fail, but
  worth a decision before the surface is frozen. Raised 2026-09-23.
- The minimum Rust version (1.85) is set by clap 4.6; CI builds with it.
- More speed in the in-memory engine. Profile for 20,000 taxa (about
  15 s): heap pushes of each new node's distances 6 s (4 s of it the
  sift-up), the update loop 4 s, pulls 2 s, rebuilds 1 s. Radix-sorted
  batches in place of the post-rebuild heaps would cut the pushes; the
  update loop is bound by cache misses on the triangular matrix.
- Reduce in-memory engine memory: the queues hold an entry for every pair
  (about 12 bytes each); 104 GB at 100,000 taxa.
- External-memory engine, 20,000 taxa at a 100 MB budget (time not yet measured on an
  unloaded machine): level
  merges of the disk heaps dominate, then pulls, the update loop, spill
  sorts and staging pushes. More slots per level would reduce how many
  merges each entry passes through.
- Branch lengths can come out negative when both children of a join get a
  negative length (the reference does the same; 93 of 40,000 at 20,000
  simulated taxa). Decide whether to clamp both to zero.
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
- The plan still under-accounts by about 6%: at 100,000 taxa it raises the budget to 121 MB
  and the run peaks at 128 MB. Raised 2026-09-24.
- [x] 2026-09-24: A 100,000-taxon run completes: 51m56s, 128 MB peak, 99,999 joins, a tree
  with 100,000 distinct labels. It had never finished before; see the freeze hang below.
- [x] 2026-09-24: Fixed the hang in the candidate-heap drain (appending from inside the
  loop over `cand_heaps` re-indexed it and, with one heap allowed, fed pairs straight back
  into the heap they came from). Regression test added; verified it hangs without the fix.
- Benchmark and probe runs must cap worker threads (`ninja -T N`, or a rayon pool in a test
  binary). The distance step otherwise takes all 192 cores on this shared machine, which
  disrupts other users. Raised 2026-09-24.
- Get a 20,000-taxon external-memory time on an unloaded machine. Runs this session ranged
  from 80 s to 445 s depending on other users' load. Raised 2026-09-22.
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
- Release 2.0.0 final once the library API review is done. Travis asked for it on
  2026-09-21 and pulled it back when he saw the open items.

- [x] 2026-09-15: `-v` with no value fails ("a value is required for --verbose"); accept a bare `-v` as verbose level 2 (clap default_missing_value). Seen by a colleague 2026-09-14.
- [x] 2026-09-16, decided against: Quote Newick labels that contain `#` (and any other character outside the plain set): extended-Newick readers such as IcyTree and Dendroscope read `name#tag` as a reticulation node and merge every leaf sharing the tag, so RepeatMasker-style names come out as a network with cycles. Found 2026-09-14. Quoting only helps IcyTree; Dendroscope and SplitsTree strip quotes and still merge labels whose tag starts with H, L or R, so Travis dropped the idea as a hack with little benefit.
- [x] 2026-09-21: `--collapse_identical` is now the default, with `--no_collapse_identical` to turn it off.
- [x] 2026-09-22: `--memory` now bounds the external-memory engine (20,000 taxa: 309 MB
  down to 76 MB). Written and tested, not yet committed.
