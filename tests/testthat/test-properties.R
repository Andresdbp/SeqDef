# Properties documented in the JEB revision (2026): exactness of the O(n)
# traversals beyond the cases in test-SeqDef.R, the Brownian/linear identity,
# clade-level ties, the returned unnormalised availability, and edge cases.

dense_A <- function(tree, s, lambda, kernel) {
  D  <- ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label]
  td <- max(ape::node.depth.edgelength(tree))
  x  <- D / td
  W <- switch(kernel,
              exponential = exp(-lambda * x),
              brownian    = {
                depth <- ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)]
                C <- outer(depth, depth, "+") / 2 - D / 2   # shared root-to-MRCA path length
                C / sqrt(outer(depth, depth))
              })
  as.numeric(W %*% s)
}

test_that("O(n) exponential traversal is exact on ultrametric trees, continuous S and polytomies", {
  set.seed(21)
  for (k in 1:5) {
    tr <- ape::rcoal(80)
    if (k > 3) tr <- ape::di2multi(tr, tol = stats::quantile(tr$edge.length, 0.3))
    df <- data.frame(tip = tr$tip.label, S = stats::runif(80))
    for (lam in c(0.5, 3.4, 20)) {
      r <- suppressMessages(SeqDef(tr, df, lambda = lam))
      expect_equal(unname(r$availability), dense_A(tr, df$S, lam, "exponential"), tolerance = 1e-10)
    }
  }
  expect_false(ape::is.binary(tr))
})

test_that("O(n) Brownian traversal equals the dense phylogenetic-correlation calculation", {
  set.seed(22)
  for (tr in list(ape::rcoal(60), ape::rtree(60))) {
    df <- data.frame(tip = tr$tip.label, S = stats::rbinom(60, 1, 0.3))
    df$S[1] <- 1
    r <- suppressMessages(SeqDef(tr, df, kernel = "brownian"))
    expect_equal(unname(r$availability), dense_A(tr, df$S, NA, "brownian"), tolerance = 1e-10)
    expect_true(is.na(r$lambda))
  }
})

test_that("on an ultrametric tree the Brownian kernel is exactly the linear kernel at lambda = 0.5", {
  set.seed(23)
  tr <- ape::rcoal(70)
  df <- data.frame(tip = tr$tip.label, S = stats::rbinom(70, 1, 0.3)); df$S[3] <- 1
  bm  <- suppressMessages(SeqDef(tr, df, kernel = "brownian"))
  lin <- suppressMessages(SeqDef(tr, df, kernel = "linear", lambda = 0.5))
  expect_equal(unname(bm$availability), unname(lin$availability), tolerance = 1e-10)
})

test_that("all tips of a clade without data receive identical scores (ultrametric tree)", {
  tr <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):2,((e:2,f:2):1,g:3):1);")
  df <- data.frame(tip = tr$tip.label, S = c(1, 0, 0, 0, 0, 0, 1))
  r <- suppressMessages(SeqDef(tr, df, lambda = 3))
  expect_equal(r$seqdef[["c"]], r$seqdef[["d"]])
  expect_equal(r$seqdef[["e"]], r$seqdef[["f"]])
  expect_lt(r$seqdef[["b"]], r$seqdef[["c"]])          # sister of a sequenced tip
})

test_that("$availability is the unnormalised A_f and $seqdef its min-max complement", {
  set.seed(24)
  tr <- ape::rcoal(30)
  df <- data.frame(tip = tr$tip.label, S = stats::rbinom(30, 1, 0.3)); df$S[1] <- 1
  r1 <- suppressMessages(SeqDef(tr, df, lambda = 2, scale = FALSE, invert = FALSE))
  r2 <- suppressMessages(SeqDef(tr, df, lambda = 2))
  expect_equal(r1$seqdef, r1$availability)
  expect_equal(r2$availability, r1$availability)
  A <- r2$availability
  expect_equal(r2$seqdef, 1 - (A - min(A)) / (max(A) - min(A)))
  expect_equal(range(r2$seqdef), c(0, 1))
  expect_identical(r2$kernel, "exponential")
  expect_equal(r2$tree_depth, max(ape::node.depth.edgelength(tr)))
})

test_that("no data anywhere: SeqDef = 1 for every tip, with a warning", {
  set.seed(25)
  tr <- ape::rcoal(20)
  df0 <- data.frame(tip = tr$tip.label, S = 0)
  expect_warning(r <- suppressMessages(SeqDef(tr, df0, lambda = 2)), "same phylogenetically weighted availability")
  expect_true(all(r$seqdef == 1))
  expect_true(all(r$availability == 0))
})

test_that("numeric lambda must be a single positive number", {
  tr <- ape::rcoal(10); df <- data.frame(tip = tr$tip.label, S = c(1, rep(0, 9)))
  expect_error(SeqDef(tr, df, lambda = -1), "positive")
  expect_error(SeqDef(tr, df, lambda = c(1, 2)), "positive")
})

test_that("large lambda reproduces the ranking by distance to the nearest tip with data", {
  set.seed(26)
  tr <- ape::rcoal(60)
  s <- stats::rbinom(60, 1, 0.15); s[1] <- 1
  df <- data.frame(tip = tr$tip.label, S = s)
  D <- ape::cophenetic.phylo(tr)[tr$tip.label, tr$tip.label]
  Ds <- D[, s == 1, drop = FALSE]; Ds[Ds == 0] <- Inf
  nsd <- apply(Ds, 1, min)
  u <- s == 0
  r <- suppressMessages(SeqDef(tr, df, lambda = 25))
  expect_gt(stats::cor(r$seqdef[u], nsd[u], method = "spearman"), 0.9)
})
