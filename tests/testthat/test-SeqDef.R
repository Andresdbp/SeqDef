# Brute-force dense reference implementation of the exponential-kernel
# availability, used to validate the O(n) traversal.
dense_exp_avail <- function(tree, s, lambda, td) {
  dmat <- ape::cophenetic.phylo(tree)
  dmat <- dmat[tree$tip.label, tree$tip.label]
  as.numeric(exp(-lambda * dmat / td) %*% s)
}

test_that("O(n) exponential traversal matches brute-force dense computation", {
  set.seed(1)
  tree <- ape::rtree(50)
  s <- sample(c(0, 1), 50, replace = TRUE)
  td <- max(ape::node.depth.edgelength(tree))
  for (lam in c(0.5, 1, 3, 10, 25)) {
    fast <- SeqDef:::.seqdef_exp_avail(tree, s, lam, td)
    ref  <- dense_exp_avail(tree, s, lam, td)
    expect_equal(unname(fast), ref, tolerance = 1e-8)
  }
})

test_that("SeqDef end-to-end matches dense computation with fixed lambda", {
  set.seed(2)
  tree <- ape::rtree(40)
  df <- data.frame(species = tree$tip.label,
                   availability = sample(c(0, 1), 40, replace = TRUE))
  res <- SeqDef(tree, df, lambda = 5, invert = FALSE, scale = FALSE)
  td <- max(ape::node.depth.edgelength(tree))
  ref <- dense_exp_avail(tree, df$availability[match(tree$tip.label, df$species)], 5, td)
  expect_equal(unname(res$seqdef), ref, tolerance = 1e-8)
  expect_named(res$seqdef, tree$tip.label)
})

test_that("auto_max warns and falls back on all-equal data (seqdef-01)", {
  set.seed(3)
  tree <- ape::rtree(20)
  df <- data.frame(species = tree$tip.label, availability = rep(1, 20))
  expect_warning(res <- suppressMessages(SeqDef(tree, df, lambda = "auto_max")),
                 "no variation", ignore.case = TRUE)
  expect_false(anyNA(res$seqdef))
  # all-zero (a clade with no genomes yet) crashed the old code with `if (NA)`
  df0 <- data.frame(species = tree$tip.label, availability = rep(0, 20))
  expect_warning(res0 <- suppressMessages(SeqDef(tree, df0, lambda = "auto_max")),
                 "no variation", ignore.case = TRUE)
  expect_false(anyNA(res0$seqdef))
})

test_that("all-equal data still works with a fixed numeric lambda", {
  set.seed(3)
  tree <- ape::rtree(20)
  df <- data.frame(species = tree$tip.label, availability = rep(1, 20))
  res <- suppressMessages(SeqDef(tree, df, lambda = 3))
  expect_false(anyNA(res$seqdef))                        # no NaN/NA anywhere
  expect_true(all(res$seqdef >= 0 & res$seqdef <= 1))    # valid scaled scores

  # A truly degenerate score vector (star tree: all tips equidistant) rescales
  # to 0 instead of NaN.
  star <- ape::stree(8, type = "star")
  star$edge.length <- rep(1, nrow(star$edge))
  df2 <- data.frame(species = star$tip.label, availability = rep(1, 8))
  res2 <- suppressMessages(SeqDef(star, df2, lambda = 3, invert = FALSE))
  expect_false(anyNA(res2$seqdef))
  expect_true(all(res2$seqdef == 0))
})

test_that("NA in data column errors informatively (seqdef-02)", {
  set.seed(4)
  tree <- ape::rtree(20)
  vals <- sample(c(0, 1), 20, replace = TRUE)
  vals[c(3, 7)] <- NA
  df <- data.frame(species = tree$tip.label, availability = vals)
  expect_error(SeqDef(tree, df, lambda = 3), "NA")
})

test_that("factor data column is converted via labels, not level codes (seqdef-04)", {
  set.seed(5)
  tree <- ape::rtree(20)
  vals <- sample(c(0, 1), 20, replace = TRUE)
  df_num <- data.frame(species = tree$tip.label, availability = vals)
  df_fac <- data.frame(species = tree$tip.label,
                       availability = factor(vals))  # levels "0","1" -> codes 1,2
  res_num <- suppressMessages(SeqDef(tree, df_num, lambda = 3))
  res_fac <- suppressMessages(SeqDef(tree, df_fac, lambda = 3))
  expect_equal(res_fac$seqdef, res_num$seqdef)
  expect_equal(res_fac$empdata, res_num$empdata)
})

test_that("non-numeric factor levels error rather than silently using codes", {
  set.seed(6)
  tree <- ape::rtree(10)
  df <- data.frame(species = tree$tip.label,
                   availability = factor(sample(c("yes", "no"), 10, replace = TRUE)))
  expect_error(SeqDef(tree, df, lambda = 3), "non-numeric")
})

test_that("empty tree/data intersection errors informatively (seqdef-06)", {
  set.seed(7)
  tree <- ape::rtree(10)
  df <- data.frame(species = gsub("t", "taxon ", tree$tip.label),  # mismatched labels
                   availability = sample(c(0, 1), 10, replace = TRUE))
  expect_error(SeqDef(tree, df, lambda = 3), "No taxa shared")
})

test_that("duplicated taxa warn and use the first row (seqdef-06)", {
  set.seed(8)
  tree <- ape::rtree(10)
  df <- data.frame(species = c(tree$tip.label, tree$tip.label[1]),
                   availability = c(sample(c(0, 1), 10, replace = TRUE), 99))
  expect_warning(res <- suppressMessages(SeqDef(tree, df, lambda = 3)), "Duplicated")
  expect_equal(unname(res$empdata[match(tree$tip.label[1], tree$tip.label)]),
               df$availability[1])
})

test_that("tree without branch lengths errors informatively (seqdef-05/06)", {
  tree <- ape::rtree(10)
  tree$edge.length <- NULL
  df <- data.frame(species = paste0("t", 1:10), availability = rep(c(0, 1), 5))
  expect_error(SeqDef(tree, df, lambda = 3), "branch lengths")
})

test_that("non-ultrametric trees use max root-to-tip depth (seqdef-05)", {
  set.seed(9)
  tree <- ape::rtree(20)  # rtree() is non-ultrametric
  expect_false(ape::is.ultrametric(tree))
  df <- data.frame(species = tree$tip.label,
                   availability = sample(c(0, 1), 20, replace = TRUE))
  res <- suppressMessages(SeqDef(tree, df, lambda = 5, invert = FALSE, scale = FALSE))
  td <- max(ape::node.depth.edgelength(tree))
  ref <- dense_exp_avail(tree, df$availability[match(tree$tip.label, df$species)], 5, td)
  expect_equal(unname(res$seqdef), ref, tolerance = 1e-8)
})

test_that("by_genus dense and subtree paths give identical lambda (seqdef-08)", {
  set.seed(10)
  tree <- ape::rtree(60)
  genera <- rep(paste0("Genus", 1:12), each = 5)
  tree$tip.label <- paste(genera, seq_len(60), sep = "_")
  td <- max(ape::node.depth.edgelength(tree))
  dense   <- SeqDef:::.seqdef_intra_genus_dists(tree, td, dense_cutoff = 4000)
  subtree <- SeqDef:::.seqdef_intra_genus_dists(tree, td, dense_cutoff = 0)
  expect_equal(sort(dense), sort(subtree), tolerance = 1e-10)
})

test_that("gaussian and linear kernels run and match direct dense computation", {
  set.seed(11)
  tree <- ape::rtree(25)
  df <- data.frame(species = tree$tip.label,
                   availability = sample(c(0, 1), 25, replace = TRUE))
  s <- df$availability[match(tree$tip.label, df$species)]
  td <- max(ape::node.depth.edgelength(tree))
  nd <- ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label] / td

  res_g <- suppressMessages(SeqDef(tree, df, lambda = 4, kernel = "gaussian",
                                   invert = FALSE, scale = FALSE))
  expect_equal(unname(res_g$seqdef), as.numeric(exp(-4 * nd^2) %*% s), tolerance = 1e-8)

  res_l <- suppressMessages(SeqDef(tree, df, lambda = 2, kernel = "linear",
                                   invert = FALSE, scale = FALSE))
  z <- 1 - 2 * nd
  expect_equal(unname(res_l$seqdef), as.numeric((z * (z > 0)) %*% s), tolerance = 1e-8)
})
