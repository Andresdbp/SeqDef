#' Calculate Sequencing Deficiency (SeqDef) Scores
#'
#' This function computes a sequencing deficiency score for tips on a phylogenetic tree.
#' It utilizes pairwise cophenetic distances and a data attribute associated with the tips
#' to calculate a "concordance" or deficiency metric.
#'
#' @param tree A phylogenetic tree object of class \code{phylo} (from the \code{ape} package).
#' @param df A data frame containing tip labels in the first column and the data of interest in \code{data.col}.
#' @param data.col An integer or character string indicating the column index or name in \code{df} containing the numeric data to analyze. Defaults to 2.
#' @param invert Logical. If \code{TRUE}, the score is calculated as \code{1 - mean(x)}. Defaults to \code{FALSE}.
#' @param scale Logical. If \code{TRUE}, the final scores are min-max scaled between 0 and 1. Defaults to \code{TRUE}.
#'
#' @return A list of class \code{"seqdef"} containing:
#' \item{tree}{The input phylogenetic tree.}
#' \item{seqdef}{A named vector of calculated SeqDef scores.}
#' \item{empdata}{The original empirical data column used for calculations.}
#'
#' @importFrom ape branching.times cophenetic.phylo
#' @export
#'
#' @examples
#' \dontrun{
#' library(ape)
#' # Create a random tree and data
#' tree <- rtree(10)
#' df <- data.frame(species = tree$tip.label, seqdat = sample(c(0,1), 10, replace = T))
#' 
#' # Run function
#' result <- SeqDef(tree, df, data.col = "seqdat")
#' }
SeqDef <- function(tree, df, data.col = 2, invert = FALSE, scale = TRUE){
  
  # Check if tip count matches data rows
  if(length(tree$tip.label) != nrow(df)){
    stop("Tree and Data should have the same length (number of tips == rows).")
  } 
  
  # Ensure column selection handles names or indices
  if(is.character(data.col)){
    data_values <- df[[data.col]]
  } else {
    data_values <- df[, data.col]
  }
  
  # Max branching time (tree depth)
  td <- max(ape::branching.times(tree))
  
  # Pairwise distances
  dist.matrix <- ape::cophenetic.phylo(tree)
  
  # Calculate distance proportion
  # (dist.matrix / 2) approximates distance to MRCA from tip
  dist.prop <- 1 - ((dist.matrix/2) / td)
  
  # Initialize results vector
  synscores <- vector(length = length(tree$tip.label))
  names(synscores) <- tree$tip.label
  
  # Double loop optimization note: This is O(N^2). 
  # For very large trees, vectorization would be preferred, 
  # but keeping original logic for fidelity.
  
  tip_labels <- tree$tip.label
  
  for(i in seq_along(tip_labels)){
    focal.tip <- tip_labels[i]
    
    # Extract the row of proportions for the focal tip
    # We can subset the matrix directly rather than looping j
    focal_props <- dist.prop[focal.tip, ]
    
    # We need to align the data_values with the order of focal_props columns
    # Assuming df order matches tree labels is risky, so we match indices:
    match_indices <- match(names(focal_props), df[[1]])
    current_data_values <- data_values[match_indices]
    
    # Calculate x vector
    x <- focal_props * current_data_values
    
    if(invert){
      synscores[i] <- 1 - mean(x, na.rm = TRUE)
    } else {
      synscores[i] <- mean(x, na.rm = TRUE)
    }
  }
  
  # Scale to 0-1
  if(scale){
    rng <- max(synscores) - min(synscores)
    if(rng == 0) rng <- 1 # Prevent division by zero if all scores are identical
    synscores <- (synscores - min(synscores)) / rng
  }
  
  results <- list(tree = tree, seqdef = synscores, empdata = data_values)
  class(results) <- "seqdef"
  
  return(results)
}
