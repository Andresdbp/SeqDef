#' Calculate Sequencing Deficiency (SeqDef) Scores
#'
#' Computes a sequencing deficiency score for tips on a phylogenetic tree.
#' It utilizes pairwise cophenetic distances and a data attribute to calculate 
#' a phylogenetically weighted "deficiency" metric.
#'
#' @param tree A phylogenetic tree object of class \code{phylo}.
#' @param df A data frame containing tip labels (col 1) and data (col 2 or \code{data.col}).
#' @param data.col Integer or character. Column index/name for the data. Defaults to 2.
#' @param invert Logical. If TRUE, returns (1 - Score). Defaults to TRUE.
#' @param scale Logical. If TRUE, min-max scales result to 0-1. Defaults to TRUE.
#' @param lambda optimization method ("auto_max", "by_genus") or a numeric value.
#'
#' @return A list of class \code{"seqdef"} containing the tree, scores, data, and lambda used.
#' @importFrom ape branching.times cophenetic.phylo keep.tip
#' @export
SeqDef <- function(tree, df, data.col = 2, invert = TRUE, scale = TRUE, lambda = "auto_max"){
  
  # 1. Convert & Align Data
  df <- as.data.frame(df)
  
  # Prune tree if it has tips not in data
  common_taxa <- intersect(tree$tip.label, df[,1])
  if(length(common_taxa) < length(tree$tip.label)){
    tree <- ape::keep.tip(tree, common_taxa)
  }
  
  # Align dataframe to match tree order exactly
  df <- df[match(tree$tip.label, df[,1]), ]
  
  # 2. Calculate Tree Depth & Matrix
  td <- max(ape::branching.times(tree))
  dist.matrix <- ape::cophenetic.phylo(tree)
  dist.matrix <- dist.matrix[tree$tip.label, tree$tip.label]
  
  # 3. Handle Lambda Selection
  final_lambda <- 0
  
  if (length(lambda) == 1 && lambda == "auto_max") {
    # --- MAXIMIZE VARIANCE LOGIC (With 10% Drop Tolerance) ---
    message("Optimizing Lambda: Searching for peak variance (Stopping if variance drops >10% below peak)...")
    
    s_vec <- df[[data.col]]
    norm_dists <- dist.matrix / td
    lambda_seq <- seq(1, 50, 0.1)
    
    best_lambda <- lambda_seq[1]
    max_var <- -1
    
    calc_var <- function(x) {
      w_mat <- exp(-x * norm_dists)
      raw_scores <- as.numeric(w_mat %*% s_vec)
      rng <- range(raw_scores)
      if(rng[2] - rng[1] < 1e-9) return(0)
      norm_scores <- (raw_scores - rng[1]) / (rng[2] - rng[1])
      return(var(1 - norm_scores))
    }
    
    max_var <- calc_var(best_lambda)
    
    for(lam in lambda_seq[-1]) {
      curr_var <- calc_var(lam)
      
      if(curr_var > max_var) {
        max_var <- curr_var
        best_lambda <- lam
      } else {
        drop_ratio <- (max_var - curr_var) / max_var
        if(drop_ratio > 0.10) {
          break
        }
      }
    }
    final_lambda <- best_lambda
    message(sprintf("Selected Lambda: %.1f (Peak Variance: %.5f)", final_lambda, max_var))
    
  } else if (length(lambda) == 1 && lambda == "by_genus") {
    # --- GENUS-SCALE LOGIC ---
    genera <- sapply(strsplit(tree$tip.label, "[_ ]"), `[`, 1)
    unique_genera <- unique(genera)
    intra_genus_dists <- c()
    
    for(g in unique_genera) {
      tips_in_genus <- tree$tip.label[genera == g]
      if(length(tips_in_genus) > 1) {
        sub_mat <- dist.matrix[tips_in_genus, tips_in_genus]
        d_vals <- sub_mat[lower.tri(sub_mat)]
        d_vals <- d_vals[d_vals > 0]
        intra_genus_dists <- c(intra_genus_dists, d_vals / td)
      }
    }
    
    if(length(intra_genus_dists) > 0){
      median_dist <- median(intra_genus_dists)
      final_lambda <- log(2) / median_dist
      message(sprintf("Genus-Scale Lambda: %.2f (Half-life: %.4f)", final_lambda, median_dist))
    } else {
      warning("Biological calibration failed: No genera with >1 species found. Defaulting to Lambda=3.")
      final_lambda <- 3
    }
    
  } else {
    if(!is.numeric(lambda)) stop("Lambda must be 'auto_max', 'by_genus', or numeric.")
    final_lambda <- lambda
  }
  
  # 4. Calculation
  dist.prop <- exp(-final_lambda * dist.matrix / td)
  raw_scores <- dist.prop %*% as.numeric(df[, data.col])
  synscores <- as.vector(raw_scores)
  names(synscores) <- tree$tip.label
  
  # 5. Scale & Invert
  if(scale){
    rng <- max(synscores, na.rm = TRUE) - min(synscores, na.rm = TRUE)
    if(rng == 0) synscores[] <- 0 else synscores <- (synscores - min(synscores, na.rm = TRUE)) / rng
  }
  if(invert) synscores <- 1 - synscores
  
  # 6. Return
  results <- list(tree = tree, seqdef = synscores, empdata = df[, data.col], lambda = final_lambda)
  class(results) <- "seqdef"
  return(results)
}