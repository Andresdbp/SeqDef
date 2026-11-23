# SeqDef 0.2.0

## New features

* `SeqDef()` now also returns the unnormalised phylogenetically weighted
  availability (`$availability`, A_f), the kernel used (`$kernel`) and the tree
  depth used to normalise distances (`$tree_depth`). Min-max scaled scores are
  relative to the tree analysed, so `$availability` is the quantity to compare
  across trees or data snapshots.

## Behaviour and documentation

* When every tip has the same availability (e.g. no tip has data) the min-max
  denominator is zero; scaled availability is set to 0 (SeqDef = 1 for every
  tip) as before, but a warning now says so.
* A numeric `lambda` must be a single positive number.
* Documentation now states: (i) on an ultrametric tree all tips of a clade
  without data receive the same score, so SeqDef ranks lineages and a second
  criterion (e.g. `calc_priority()`) separates species within them;
  (ii) normalised distances d/T lie in [0, 2]; (iii) the Brownian kernel equals
  the linear kernel with lambda = 0.5 on an ultrametric tree; (iv) large lambda
  reproduces the ranking by distance to the nearest tip with data; (v) min-max
  scaling is a within-tree rescaling; (vi) how `calc_priority()`'s base
  controls the weight of the trait.

## Tests

* New tests: O(n) exponential and Brownian traversals versus the full-matrix
  calculation (ultrametric and non-ultrametric trees, continuous data,
  polytomies), the Brownian/linear(0.5) identity, clade-level ties, the
  returned availability, the no-data case, and the large-lambda limit.

# SeqDef 0.1.1

## Bug fixes

* `SeqDef()` no longer crashes with an obscure `if (NA)` error when the data
  column has no variation (e.g. a clade with no genomes yet): `auto_max` now
  stops with an informative message suggesting a numeric lambda or
  `by_genus`, and the variance-drop stopping rule guards against a zero peak
  variance.
* `SeqDef()` validates its inputs up front: NA or non-numeric values in the
  data column are reported by taxon name instead of poisoning every score;
  factor data columns are converted via their labels rather than their
  internal level codes; trees without branch lengths, and tree/data pairs
  with no shared taxa (e.g. spaces-vs-underscores label mismatches), error
  informatively; duplicated taxa in `df` warn and use the first row; silently
  dropped tips are now reported.
* Tree depth is now computed with `max(node.depth.edgelength(tree))` instead
  of `max(branching.times(tree))`, which is only correct for ultrametric
  trees. For non-ultrametric trees this changes the lambda scale (intended).
* `calc_priority()` now errors on unnamed trait vectors (previously it
  silently returned scores identical to the SeqDef scores) and warns with a
  count and examples when taxa are missing trait values.
* Degeneracy tolerance is now consistently `< 1e-9` in both the lambda
  search and the final rescaling; `$empdata` now stores the numeric vector
  actually used for scoring (in tip order) rather than the raw column.

## Performance

* `by_genus` lambda calibration is now size-aware: trees up to 4000 tips use
  a single dense cophenetic matrix (fastest at small n), while larger trees
  compute intra-genus distances from per-genus subtrees, avoiding the two
  full n-by-n matrix allocations that exhausted memory at 10k+ tips. A
  redundant second full-matrix copy was removed from the dense path.
* Lambda-independent preparation (postorder edge reordering) is hoisted out
  of the `auto_max` per-lambda loop for the exponential kernel, and the
  Gaussian kernel's squared-distance matrix is computed once instead of per
  lambda iteration.

## Infrastructure

* Added the missing LICENSE file (MIT).
* Added a testthat suite (kernel equivalence against a brute-force dense
  implementation, degenerate-data guards, `calc_priority` validation, factor
  handling, by_genus dense/subtree equivalence).
* DESCRIPTION cleanup: removed `LazyData` (no data sets), removed
  `VignetteBuilder` (no vignettes), replaced the `tidyverse` meta-package in
  Suggests with the packages actually used; `stats` functions (`var`,
  `median`, `setNames`) are now properly imported.
