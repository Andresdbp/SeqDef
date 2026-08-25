# Linear-time, linear-memory availability for the EXPONENTIAL kernel.
#
# A_f = sum_i exp(-lambda * d_fi / td) * s_i, where d_fi is the cophenetic (path)
# distance. Because d_fi is additive along the tree path and the kernel is
# exponential, exp(-lambda d_fi/td) factorizes into a product of per-edge factors
# phi_e = exp(-lambda L_e / td). A_f is then computed by a two-pass sum-product
# traversal (Felsenstein-style message passing) in O(n) time and O(n) memory --
# no n x n distance matrix is formed. This is EXACT (identical to the dense
# computation) and is specific to the exponential kernel (the only one of the
# three that factorizes into independent per-branch factors).

# Precompute the lambda-independent parts of the exponential traversal
# (postorder edge ordering). Hoisting this out of the per-lambda loop avoids
# repeating reorder.phylo/name matching on every lambda evaluation.
.seqdef_exp_prep <- function(tree) {
  n   <- length(tree$tip.label)
  phy <- ape::reorder.phylo(tree, "postorder")   # children before parents
  list(n = n, E = phy$edge, elen = phy$edge.length,
       nN = n + phy$Nnode, tip.label = phy$tip.label)
}

.seqdef_exp_avail_prep <- function(prep, s, lambda, td) {
  n   <- prep$n
  E   <- prep$E
  phi <- exp(-lambda * prep$elen / td)            # per-edge transmission factor

  # Up-sweep (post-order): Down[v] = sum over leaves below v of s_i * (product of
  # phi from v down to i). Leaves: Down = s. Internal: sum_children phi * Down.
  Down <- numeric(prep$nN)
  Down[seq_len(n)] <- s
  for (k in seq_len(nrow(E))) {
    Down[E[k, 1]] <- Down[E[k, 1]] + phi[k] * Down[E[k, 2]]
  }

  # Down-sweep (pre-order): Up[v] = signal from everything NOT below v, at v.
  # Up[child] = phi * (Up[parent] + Down[parent] - phi * Down[child]).
  Up <- numeric(prep$nN)                          # Up[root] = 0
  for (k in rev(seq_len(nrow(E)))) {
    p <- E[k, 1]; c <- E[k, 2]
    Up[c] <- phi[k] * (Up[p] + Down[p] - phi[k] * Down[c])
  }

  A <- Down[seq_len(n)] + Up[seq_len(n)]          # tip f: A_f = s_f + Up[f]
  names(A) <- prep$tip.label                      # tip numbering unchanged by reorder
  A
}

.seqdef_exp_avail <- function(tree, s, lambda, td) {
  .seqdef_exp_avail_prep(.seqdef_exp_prep(tree), s, lambda, td)
}

# Brownian-motion covariance kernel (parameter-free), correlation form, O(n).
#
# Weight w_fi = phylogenetic correlation under BM = C_fi / sqrt(C_ff C_ii), with
# C_fi = shared root-to-MRCA path length and C_ii = root-to-tip depth (w in
# [0,1], diagonal 1). There is no lambda. Availability A_f = sum_i w_fi s_i is
# computed in O(n) via the three-point structure: with s'_i = s_i / sqrt(C_ii)
# and subtree sums S'(v), A_f = (1/sqrt(C_ff)) sum_{v on root->f path} L_v S'(v).
# For an ultrametric tree this equals the exponential kernel's flat limit,
# w_fi = 1 - d_fi/(2T) (i.e. the exponential at a fixed small lambda ~ 0.5).
.seqdef_bm_avail <- function(tree, s) {
  n   <- length(tree$tip.label)
  phy <- ape::reorder.phylo(tree, "postorder")
  E   <- phy$edge; L <- phy$edge.length; nN <- n + phy$Nnode
  depth <- ape::node.depth.edgelength(phy)         # root-to-node distance
  dtip  <- depth[seq_len(n)]; dtip[dtip <= 0] <- .Machine$double.eps
  sp    <- s / sqrt(dtip)                          # s'_i = s_i / sqrt(C_ii)

  Sp <- numeric(nN); Sp[seq_len(n)] <- sp          # subtree sums of s'
  for (k in seq_len(nrow(E))) Sp[E[k, 1]] <- Sp[E[k, 1]] + Sp[E[k, 2]]

  G <- numeric(nN)                                 # G[v] = sum_{u on root->v} L_u S'(u)
  for (k in rev(seq_len(nrow(E)))) {
    p <- E[k, 1]; c <- E[k, 2]
    G[c] <- G[p] + L[k] * Sp[c]
  }
  A <- G[seq_len(n)] / sqrt(dtip)
  names(A) <- phy$tip.label
  A[tree$tip.label]
}

# Intra-genus normalized cophenetic distances for the "by_genus" calibration.
# Size-aware: up to `dense_cutoff` tips (or when a dense cophenetic matrix is
# already available from the kernel computation) a single dense matrix is
# fastest; above the cutoff the O(n^2) memory is prohibitive, so distances are
# computed per genus on extracted subtrees (identical values, O(n) memory).
.seqdef_intra_genus_dists <- function(tree, td, dist.matrix = NULL,
                                      dense_cutoff = 4000) {
  genera <- sapply(strsplit(tree$tip.label, "[_ ]"), `[`, 1)
  use_dense <- !is.null(dist.matrix) ||
    length(tree$tip.label) <= dense_cutoff
  if (use_dense && is.null(dist.matrix)) {
    dist.matrix <- ape::cophenetic.phylo(tree)
  }
  intra_genus_dists <- c()
  for (g in unique(genera)) {
    tips_in_genus <- tree$tip.label[genera == g]
    if (length(tips_in_genus) > 1) {
      if (use_dense) {
        sub_mat <- dist.matrix[tips_in_genus, tips_in_genus]
      } else {
        sub_mat <- ape::cophenetic.phylo(ape::keep.tip(tree, tips_in_genus))
      }
      d_vals <- sub_mat[lower.tri(sub_mat)]
      d_vals <- d_vals[d_vals > 0]
      intra_genus_dists <- c(intra_genus_dists, d_vals / td)
    }
  }
  intra_genus_dists
}

#' Calculate Sequencing Deficiency (SeqDef) Scores
#'
#' Computes a sequencing deficiency score for tips on a phylogenetic tree.
#' It utilizes pairwise cophenetic distances and a data attribute to calculate
#' a phylogenetically weighted "deficiency" metric.
#'
#' For the (default) exponential kernel the phylogenetically weighted
#' availability is computed by an exact O(n)-time, O(n)-memory tree traversal
#' (the exponential of an additive distance factorizes into per-branch factors),
#' so no n-by-n distance matrix is formed and the method scales to very large
#' trees. The Gaussian and linear kernels use the dense O(n^2) cophenetic matrix.
#' The Brownian-motion kernel is a parameter-free, O(n) alternative based on the
#' phylogenetic correlation (shared ancestry); it is the flat, low-lambda limit
#' of the exponential kernel.
#'
#' Distances are normalized by the tree depth (the maximum root-to-tip path
#' length), which is well defined for both ultrametric and non-ultrametric
#' trees. Input data are validated: the tree must have branch lengths, the
#' data column must be numeric (factors are converted via their labels, never
#' their level codes) and free of NAs, taxon labels in \code{df} must overlap
#' \code{tree$tip.label}, and duplicated taxa or dropped tips are reported.
#'
#' @param tree A phylogenetic tree object of class \code{phylo} with branch lengths.
#' @param df A data frame containing tip labels (col 1) and data (col 2 or \code{data.col}).
#' @param data.col Integer or character. Column index/name for the data. Defaults to 2.
#' @param invert Logical. If TRUE, returns (1 - Score). Defaults to TRUE.
#' @param scale Logical. If TRUE, min-max scales result to 0-1. Defaults to TRUE.
#' @param lambda optimization method ("auto_max", "by_genus") or a numeric value.
#' @param kernel Weighting kernel: \code{"exponential"} (default; distance-decay,
#'   exact O(n) traversal), \code{"gaussian"} or \code{"linear"} (distance-decay,
#'   dense O(n^2)), or \code{"brownian"} (parameter-free Brownian-motion
#'   phylogenetic correlation, O(n); \code{lambda} is ignored). The default
#'   reproduces the original behaviour.
#'
#' @return A list of class \code{"seqdef"} containing the tree, scores, data
#'   (the numeric vector actually used for scoring, in tip order), and lambda used.
#' @importFrom ape cophenetic.phylo keep.tip reorder.phylo node.depth.edgelength
#' @importFrom stats var median
#' @export
SeqDef <- function(tree, df, data.col = 2, invert = TRUE, scale = TRUE, lambda = "auto_max",
                   kernel = c("exponential", "gaussian", "linear", "brownian")){

  kernel <- match.arg(kernel)

  # 1. Validate & Align Data
  if (!inherits(tree, "phylo")) stop("'tree' must be an object of class 'phylo'.")
  if (is.null(tree$edge.length)) {
    stop("'tree' has no branch lengths; SeqDef requires a tree with edge lengths.")
  }
  df <- as.data.frame(df)
  if (anyDuplicated(df[, 1]) > 0) {
    warning("Duplicated taxon labels in df[, 1]; only the first row for each taxon is used.")
    df <- df[!duplicated(df[, 1]), ]
  }
  common_taxa <- intersect(tree$tip.label, df[, 1])
  if (length(common_taxa) == 0) {
    stop("No taxa shared between tree$tip.label and df[, 1]. ",
         "Check that labels use the same formatting (e.g. spaces vs underscores).")
  }
  if (length(common_taxa) < length(tree$tip.label)) {
    message(sprintf("Dropping %d tip(s) not present in df[, 1].",
                    length(tree$tip.label) - length(common_taxa)))
    tree <- ape::keep.tip(tree, common_taxa)
  }
  df <- df[match(tree$tip.label, df[, 1]), ]
  s_raw <- df[[data.col]]
  if (is.factor(s_raw)) s_raw <- as.character(s_raw)  # use labels, never level codes
  s_vec <- suppressWarnings(as.numeric(s_raw))
  if (anyNA(s_vec)) {
    bad <- tree$tip.label[is.na(s_vec)]
    stop(sprintf("The data column contains %d NA or non-numeric value(s) (e.g. %s); SeqDef requires complete numeric data.",
                 length(bad), paste(bad[seq_len(min(5, length(bad)))], collapse = ", ")))
  }

  # 2. Tree depth (normalizes distances to a relative scale). The maximum
  #    root-to-tip depth is correct for both ultrametric and non-ultrametric
  #    trees (branching.times() assumes ultrametric trees).
  td <- max(ape::node.depth.edgelength(tree))
  if (!is.finite(td) || td <= 0) {
    stop("Tree depth is zero or undefined; check the tree's edge lengths.")
  }

  # 3. Dense matrices only when the kernel needs them (gaussian/linear).
  #    The exponential and brownian paths are matrix-free (O(n)).
  is_bygenus  <- (length(lambda) == 1 && is.character(lambda) && lambda == "by_genus")
  need_matrix <- kernel %in% c("gaussian", "linear")
  dist.matrix   <- NULL
  norm_dists    <- NULL
  norm_dists_sq <- NULL
  if (need_matrix) {
    dist.matrix <- ape::cophenetic.phylo(tree)     # already in tip.label order
    norm_dists  <- dist.matrix / td
    if (kernel == "gaussian") norm_dists_sq <- norm_dists^2  # hoisted out of the lambda loop
  }
  exp_prep <- if (kernel == "exponential") .seqdef_exp_prep(tree) else NULL

  # Raw phylogenetically-weighted availability A for a given lambda.
  avail <- function(lam) {
    if (kernel == "exponential") {
      .seqdef_exp_avail_prep(exp_prep, s_vec, lam, td)       # O(n) traversal
    } else if (kernel == "brownian") {
      .seqdef_bm_avail(tree, s_vec)                          # O(n), parameter-free
    } else if (kernel == "gaussian") {
      as.numeric(exp(-lam * norm_dists_sq) %*% s_vec)        # O(n^2) dense
    } else {                                                 # linear, clamped at 0
      z <- 1 - lam * norm_dists
      as.numeric((z * (z > 0)) %*% s_vec)
    }
  }

  # 4. Lambda selection
  final_lambda <- 0
  if (kernel == "brownian") {
    if (is.numeric(lambda)) message("Brownian kernel is parameter-free; the supplied lambda is ignored.")
    final_lambda <- NA_real_
  } else if (length(lambda) == 1 && is.character(lambda) && lambda == "auto_max") {
    lambda_seq <- seq(1, 50, 0.1)
    if (length(unique(s_vec)) < 2) {
      warning("The data column has no variation (all values identical), so lambda ",
              "cannot be optimized by variance ('auto_max'); using lambda = ",
              lambda_seq[1], ". Supply a numeric lambda or use lambda = 'by_genus' ",
              "to control this directly.")
      final_lambda <- lambda_seq[1]
    } else {
    message("Optimizing Lambda: Searching for peak variance (Stopping if variance drops >10% below peak)...")

    calc_var <- function(x) {
      raw <- avail(x)
      rng <- range(raw)
      if (rng[2] - rng[1] < 1e-9) return(0)
      var(1 - (raw - rng[1]) / (rng[2] - rng[1]))
    }

    best_lambda <- lambda_seq[1]
    max_var <- calc_var(best_lambda)
    for (lam in lambda_seq[-1]) {
      curr_var <- calc_var(lam)
      if (curr_var > max_var) {
        max_var <- curr_var
        best_lambda <- lam
      } else if (max_var > 0 && (max_var - curr_var) / max_var > 0.10) {
        break
      }
    }
    final_lambda <- best_lambda
    message(sprintf("Selected Lambda: %.1f (Peak Variance: %.5f)", final_lambda, max_var))
    }

  } else if (is_bygenus) {
    intra_genus_dists <- .seqdef_intra_genus_dists(tree, td, dist.matrix = dist.matrix)
    if (length(intra_genus_dists) > 0) {
      median_dist <- median(intra_genus_dists)
      final_lambda <- log(2) / median_dist
      message(sprintf("Genus-Scale Lambda: %.2f (Half-life: %.4f)", final_lambda, median_dist))
    } else {
      warning("Biological calibration failed: No genera with >1 species found. Defaulting to Lambda=3.")
      final_lambda <- 3
    }

  } else {
    if (!is.numeric(lambda)) stop("Lambda must be 'auto_max', 'by_genus', or numeric.")
    final_lambda <- lambda
  }

  # 5. Final availability + scale + invert
  synscores <- avail(final_lambda)
  names(synscores) <- tree$tip.label
  if (scale) {
    rng <- max(synscores, na.rm = TRUE) - min(synscores, na.rm = TRUE)
    if (rng < 1e-9) synscores[] <- 0 else synscores <- (synscores - min(synscores, na.rm = TRUE)) / rng
  }
  if (invert) synscores <- 1 - synscores

  results <- list(tree = tree, seqdef = synscores, empdata = s_vec, lambda = final_lambda)
  class(results) <- "seqdef"
  return(results)
}
