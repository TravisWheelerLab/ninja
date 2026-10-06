# Changelog

## 2.1.0 (2026-10-06)

* ninja clamps a negative branch length to zero and reports how many it
  clamped. Neighbor joining produces one when the distance between the pair
  being joined is itself negative, which the update formula allows on data
  that is not additive. The engines already zeroed a negative branch and
  moved its length onto the sibling, but when the pair distance was negative
  there was nowhere left to move it, and the sibling stayed negative. A
  negative edge is not a representable branch, and a downstream step that
  sums path lengths will read it and get a shorter path than the tree
  implies. On the 700-taxon test alignment this
  changes six internal branches, between -0.0005 and -0.011, to zero.
* Every public item has a documented example: 97 new doctests, covering all
  91 public items that lacked one. The two engines, the readers and writers,
  the distance calculators and the substitution matrices each show a small
  working call with its real output.

## 2.0.0 (2026-10-05)

Since rc.4; the entries below cover the rest of 2.0.

* ninja counts memory the plan was leaving out: the per-heap staging
  buffers, the join loop's per-taxon scratch, and the parts of a candidate
  slot outside the list itself. At 100,000 taxa under `--memory 0.1` it
  used to raise the budget to 121 MB and then peak at 128 MB; it now plans
  124 MB and peaks at 126.1 MB. The remaining 1.7% looks like allocator
  slack and the leaf-name strings, which the plan does not control.
* ninja reports an unknown `--matrix` name. `--matrix BLOSUM80` used to
  report that the file BLOSUM80 does not exist. The message named neither
  of the two matrices that would have worked. ninja now reads a spec with
  no separator and no extension as a name, and reports an unknown one with
  the built-ins listed.
* The external-memory engine is internal. `nj::extmem` is no longer public,
  nor are the `DiskMatrix` and `MemoryPlan` types its signature needed.
  `run` moves to that engine on its own when the matrix will not fit the
  budget, which is how the binary reaches it. Leaving `MemoryPlan`'s fields
  public would have made any later change to them a breaking change.

## 2.0.0-rc.4 (2026-09-25)

* The `--memory` budget now bounds what the external-memory engine uses.
  The engine used to size three structures without multiplying by how many
  of each it allocated, so a 100 MB budget took 309 MB at 20,000 taxa and
  1.06 GB at 100,000. At 20,000 taxa the budget now holds, at 76 MB. At
  100,000 taxa 100 MB is not achievable: the tree, the index structures and
  the narrowest usable window need 121 MB between them, so ninja raises the
  budget to that, and the run peaks at 128 MB. A budget below 100 MB it
  raises to 100 MB, with a warning.
* Fixed a hang that stopped large runs finishing. The engine drained the
  frozen candidate heaps and appended each surviving pair back to the
  candidate list, and appending can freeze that list, which merges the
  oldest heaps away and pushes a new one. Doing that from inside the loop
  over those heaps re-indexed the vector underneath it, and where the budget
  allowed a single heap the pair just taken out went straight back into the
  heap it came from. A 100,000-taxon run never finished, spinning on one
  core without completing another join; it now takes 52 minutes. ninja
  collects the candidates and appends them once the loop has finished. Only
  a candidate list long enough to freeze reaches this code, which takes
  around 100,000 taxa at a small budget, so no test covered it; one does
  now.
* ninja raises a budget too small for the taxon count. It used to warn that
  it needed more and carry on regardless, so `--memory` bounded nothing at
  the point where a bound was worth having. It now raises the budget to what
  the tree, the index structures and the narrowest usable window need, and
  keeps the cluster count it had settled on, since a larger budget would
  otherwise buy more cluster-pair heaps and pull the requirement up again.
* Budget warnings print even under `--quiet`. The test meant to keep them
  visible fired only when the budget had been raised to the 100 MB minimum.
  A run that needed more memory than it was given therefore said nothing,
  which is the one case where the warning matters.
* A tight budget buys fewer clusters: rather than starve every cluster-pair
  heap, the engine reduces the cluster count until each heap has a megabyte
  to work with, and says so under `--verbose`.
* The library's public surface is narrower. The priority-queue module, the
  gap-counting helpers, the disk-matrix accessors and layout constants, the
  candidate heap and `cluster::write_table` are now internal; most were
  reachable only because the module holding them was left public. `NjStats`
  joins `Method` and `NjParams` at the crate root: it is what
  `nj::inmem::build` returns, and it was the only one of the three a caller
  had to name by its full path.
* Four functions turned out to be uncalled and are gone:
  `DiskMatrix::prefetch`, `DiskMatrix::mem_row`, `MinHeap::iter` and
  `CandidateHeap::is_empty`. Deleting the first removed the crate's second
  `unsafe` block, leaving the single cache-prefetch hint in
  `distance::matrix` that `lib.rs` describes.

* Identical sequences are collapsed by default: the tree is built over one
  representative of each set and the duplicates hang off it as zero-length
  branches. `--no_collapse_identical` restores the old behaviour. Only the
  arrangement of zero-length branches changes, and branch lengths shift
  slightly because the neighbor-joining formula depends on the taxon count.

## 2.0.0-rc.3 (2026-09-16)

* Protein distances use BLOSUM62 by default. `--matrix` selects BLOSUM45
  or reads a matrix file in the NCBI format; the dissimilarities are derived
  with the transformation FastTree applied to BLOSUM45, which reproduces
  its table. Protein trees therefore differ from rc.2 unless
  `--matrix BLOSUM45` is given.
* A bare `-v` now works: `-v` is verbosity 2 and `-vv` is 3. `--verbose N`
  still sets a level directly.

## 2.0.0-rc.2 (2026-09-14)

* ninja warns when sequence names contain `#`, because viewers that read
  extended Newick take `name#tag` as a reticulation node and draw the tree
  as a network with cycles.

## 2.0.0-rc.1 (2026-09-14)

NINJA's earlier releases were the 1.x line; the Rust tool starts the 2.x
line, and this is the first candidate for 2.0.0. The command line and the
library API may still change before the final release. Version 0.1.0 was
the first upload of this code to crates.io and is yanked.

* Records that share a name no longer produce a Newick file with repeated
  labels: ninja renames later occurrences `name_2`, `name_3`, ... and lists
  them on standard error. `--duplicate_names error` refuses the input instead.
* The Phylip writer formats rows in parallel and no longer caches the whole
  matrix first (31 s and 2.3 GB down to 4 s and 0.8 GB at 20,000 taxa).
  Rows no longer end with a space, which fastphylo's reader rejected.
* `--method auto` estimates the in-memory engine's size at 22 bytes per
  pair instead of 7, so a 2 GB budget now selects the external engine at
  20,000 taxa.

## 0.1.0 (2026-09-09)

First Rust release, ported from NINJA 1.2.2 (Java) with the C++ port as a
second reference.

* In-memory and external-memory neighbor-joining engines producing the same
  trees as the Java implementation on all test inputs.
* DNA and protein distances from FASTA alignments, with Jukes-Cantor,
  Kimura two-parameter, scoredist, or no correction; distance computation
  is parallel.
* Phylip distance matrices as input or output.
* From the C++ `cluster` branch: single-linkage clustering
  (`--out_type c`, `--cluster_cutoff`), Mothur's onegap distance
  (`--corr_type m`), and `--collapse_identical`, which is planned to
  become the default.
* In-memory engine about twice as fast as a direct port: sorted runs
  instead of heaps for rebuilt entries, parallel rebuilds, prefetching in
  the update loop. `--reference_order` reproduces the Java tie order.
* External-memory engine: rebuilds write sorted runs to disk directly,
  spills use the standard sort, level merges use a k-way heap, and the
  disk matrix is filled in one pass. 20,000 taxa at a 50 MB budget: 114 s
  and 1.7 GB before, 68 s and 0.3 GB after. Row sums and the criterion are
  double precision, which brings this engine's trees into agreement with
  the exact in-memory engine on every split at 20,000 taxa.
* `--method auto` names the engine choice that the Java and C++ tools
  called `default`; `default` is still accepted.
* Library API alongside the `ninja` binary.
* Integration tests against stored Java outputs; see `docs/testing.md`.
