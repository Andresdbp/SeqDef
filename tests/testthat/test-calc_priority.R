make_seqdef_res <- function(n = 10) {
  set.seed(42)
  tree <- ape::rtree(n)
  df <- data.frame(species = tree$tip.label,
                   availability = sample(c(0, 1), n, replace = TRUE))
  suppressMessages(SeqDef(tree, df, lambda = 3))
}

test_that("calc_priority errors on unnamed trait vectors (seqdef-03)", {
  res <- make_seqdef_res()
  expect_error(calc_priority(res, rep(1, length(res$seqdef))), "named")
})

test_that("calc_priority warns on missing taxa and fills with na.fill (seqdef-03)", {
  res <- make_seqdef_res()
  taxa <- names(res$seqdef)
  traits <- setNames(rep(2, length(taxa) - 3), taxa[seq_len(length(taxa) - 3)])
  expect_warning(p <- calc_priority(res, traits), "no trait value")
  # missing taxa filled with na.fill = 0 -> exponential model gives base^0 = 1
  missing <- setdiff(taxa, names(traits))
  expect_equal(p[missing], res$seqdef[missing])
  expect_equal(p[names(traits)], res$seqdef[names(traits)] * 2^2)
})

test_that("calc_priority accepts a data frame of traits", {
  res <- make_seqdef_res()
  taxa <- names(res$seqdef)
  tdf <- data.frame(species = taxa, weight = seq_along(taxa))
  p <- calc_priority(res, tdf, model = "linear")
  expect_equal(unname(p), unname(res$seqdef * seq_along(taxa)))
  expect_named(p, taxa)
})

test_that("calc_priority models behave as documented", {
  res <- make_seqdef_res()
  taxa <- names(res$seqdef)
  traits <- setNames(rep(1, length(taxa)), taxa)
  expect_equal(calc_priority(res, traits, model = "exponential", base = 2),
               res$seqdef * 2)
  expect_equal(calc_priority(res, traits, model = "additive"),
               res$seqdef + 1)
})
