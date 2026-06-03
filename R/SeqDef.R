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
.seqdef_exp_avail <- function(tree, s, lambda, td) {
  n   <- length(tree$tip.label)
  phy <- ape::reorder.phylo(tree, "postorder")   # children before parents
  E   <- phy$edge
  phi <- exp(-lambda * phy$edge.length / td)      # per-edge transmission factor
  nN  <- n + phy$Nnode

  # Up-sweep (post-order): Down[v] = sum over leaves below v of s_i * (product of
  # phi from v down to i). Leaves: Down = s. Internal: sum_children phi * Down.
  Down <- numeric(nN)
  Down[seq_len(n)] <- s
  for (k in seq_len(nrow(E))) {
    Down[E[k, 1]] <- Down[E[k, 1]] + phi[k] * Down[E[k, 2]]
  }

  # Down-sweep (pre-order): Up[v] = signal from everything NOT below v, at v.
  # Up[child] = phi * (Up[parent] + Down[parent] - phi * Down[child]).
  Up <- numeric(nN)                                # Up[root] = 0
  for (k in rev(seq_len(nrow(E)))) {
    p <- E[k, 1]; c <- E[k, 2]
    Up[c] <- phi[k] * (Up[p] + Down[p] - phi[k] * Down[c])
  }

  A <- Down[seq_len(n)] + Up[seq_len(n)]           # tip f: A_f = s_f + Up[f]
  names(A) <- phy$tip.label
  A[tree$tip.label]                                # input tip order
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
#'
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param df A data frame containing tip labels (col 1) and data (col 2 or \code{data.col}).
#' @param data.col Integer or character. Column index/name for the data. Defaults to 2.
#' @param invert Logical. If TRUE, returns (1 - Score). Defaults to TRUE.
#' @param scale Logical. If TRUE, min-max scales result to 0-1. Defaults to TRUE.
#' @param lambda optimization method ("auto_max", "by_genus") or a numeric value.
#' @param kernel Distance-decay kernel applied to the normalized cophenetic
#'   distance: \code{"exponential"} (default; uses the fast O(n) traversal),
#'   \code{"gaussian"}, or \code{"linear"}. The default reproduces the original
#'   behaviour.
#'
#' @return A list of class \code{"seqdef"} containing the tree, scores, data, and lambda used.
#' @importFrom ape branching.times cophenetic.phylo keep.tip reorder.phylo
#' @export
SeqDef <- function(tree, df, data.col = 2, invert = TRUE, scale = TRUE, lambda = "auto_max",
                   kernel = c("exponential", "gaussian", "linear")){

  kernel <- match.arg(kernel)
  # Distance-decay kernel on normalized distance x = d / tree_depth (dense path).
  kern <- function(x, lam, type) {
    switch(type,
           exponential = exp(-lam * x),
           gaussian    = exp(-lam * x^2),
           linear      = { z <- 1 - lam * x; z * (z > 0) })  # clamp at 0, preserve matrix dims
  }

  # 1. Convert & Align Data
  df <- as.data.frame(df)
  common_taxa <- intersect(tree$tip.label, df[, 1])
  if (length(common_taxa) < length(tree$tip.label)) {
    tree <- ape::keep.tip(tree, common_taxa)
  }
  df <- df[match(tree$tip.label, df[, 1]), ]
  s_vec <- as.numeric(df[[data.col]])

  # 2. Tree depth (normalizes distances to a relative scale)
  td <- max(ape::branching.times(tree))

  # 3. Distance matrix only when needed: non-exponential kernels OR by_genus.
  #    The exponential path is matrix-free (O(n)).
  is_bygenus  <- (length(lambda) == 1 && is.character(lambda) && lambda == "by_genus")
  need_matrix <- (kernel != "exponential") || is_bygenus
  norm_dists  <- NULL
  if (need_matrix) {
    dist.matrix <- ape::cophenetic.phylo(tree)
    dist.matrix <- dist.matrix[tree$tip.label, tree$tip.label]
    norm_dists  <- dist.matrix / td
  }

  # Raw phylogenetically-weighted availability A for a given lambda.
  avail <- function(lam) {
    if (kernel == "exponential")
      .seqdef_exp_avail(tree, s_vec, lam, td)                 # O(n) traversal
    else
      as.numeric(kern(norm_dists, lam, kernel) %*% s_vec)     # O(n^2) dense
  }

  # 4. Lambda selection
  final_lambda <- 0
  if (length(lambda) == 1 && is.character(lambda) && lambda == "auto_max") {
    message("Optimizing Lambda: Searching for peak variance (Stopping if variance drops >10% below peak)...")
    lambda_seq <- seq(1, 50, 0.1)

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
      } else if ((max_var - curr_var) / max_var > 0.10) {
        break
      }
    }
    final_lambda <- best_lambda
    message(sprintf("Selected Lambda: %.1f (Peak Variance: %.5f)", final_lambda, max_var))

  } else if (is_bygenus) {
    genera <- sapply(strsplit(tree$tip.label, "[_ ]"), `[`, 1)
    intra_genus_dists <- c()
    for (g in unique(genera)) {
      tips_in_genus <- tree$tip.label[genera == g]
      if (length(tips_in_genus) > 1) {
        sub_mat <- dist.matrix[tips_in_genus, tips_in_genus]
        d_vals <- sub_mat[lower.tri(sub_mat)]
        d_vals <- d_vals[d_vals > 0]
        intra_genus_dists <- c(intra_genus_dists, d_vals / td)
      }
    }
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
    if (rng == 0) synscores[] <- 0 else synscores <- (synscores - min(synscores, na.rm = TRUE)) / rng
  }
  if (invert) synscores <- 1 - synscores

  results <- list(tree = tree, seqdef = synscores, empdata = df[, data.col], lambda = final_lambda)
  class(results) <- "seqdef"
  return(results)
}
