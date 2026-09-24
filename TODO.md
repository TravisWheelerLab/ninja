# To do

Planned work, roughly in order.

- **Library use.** Make the crate a first-class library, not just the
  binary's internals: settle the public API (`DistanceCalculator`,
  `DistanceMatrix`, the two `build` functions, `Tree`, `run`), keep it
  stable across releases, document every public item with examples,
  and decide whether to offer C and Python bindings. The current
  `pub` surface exists but has not been reviewed with outside callers in
  mind.
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
- Measure what reducing the cluster count under a tight budget costs. It drops from 30 to 7
  at 20,000 taxa and to 4 at 100,000; fewer clusters mean looser bounds and more work per
  join, and nobody has measured it. Raised 2026-09-22.
- Get a 20,000-taxon external-memory time on an unloaded machine. Runs this session ranged
  from 80 s to 445 s depending on other users' load. Raised 2026-09-22.
- Find a way to run multi-hour jobs that the harness watchdog will not kill. It killed the
  100,000-taxon run twice for "low memory" while the machine had 572 GB available. A
  `setsid nohup` relaunch on 2026-09-23 is the current attempt; tmux and systemd-run are
  also on the machine and untried. Raised 2026-09-23.
- Release 2.0.0 final once the library API review is done. Travis asked for it on
  2026-09-21 and pulled it back when he saw the open items.

- [x] 2026-09-15: `-v` with no value fails ("a value is required for --verbose"); accept a bare `-v` as verbose level 2 (clap default_missing_value). Seen by a colleague 2026-09-14.
- [x] 2026-09-16, decided against: Quote Newick labels that contain `#` (and any other character outside the plain set): extended-Newick readers such as IcyTree and Dendroscope read `name#tag` as a reticulation node and merge every leaf sharing the tag, so RepeatMasker-style names come out as a network with cycles. Found 2026-09-14. Quoting only helps IcyTree; Dendroscope and SplitsTree strip quotes and still merge labels whose tag starts with H, L or R, so Travis dropped the idea as a hack with little benefit.
- [x] 2026-09-21: `--collapse_identical` is now the default, with `--no_collapse_identical` to turn it off.
- [x] 2026-09-22: `--memory` now bounds the external-memory engine (20,000 taxa: 309 MB
  down to 76 MB). Written and tested, not yet committed.
