#' Calculate Composite Priority Scores
#'
#' Combines SeqDef scores with external weighting variables (e.g., IUCN status).
#'
#' @param seqdef_res A list containing \code{$seqdef} scores (output from \code{SeqDef}).
#' @param trait_values Named vector (names matching taxa in \code{seqdef_res$seqdef})
#'   or a data frame with taxon names in column 1 and weights in column 2.
#' @param model "exponential" (default), "linear", or "additive".
#' @param base Base for exponential model (default 2).
#' @param na.fill Value to replace NAs with (default 0).
#'
#' @return A named numeric vector of priority scores.
#' @importFrom stats setNames
#' @export
calc_priority <- function(seqdef_res, trait_values, model = c("exponential", "linear", "additive"), 
                          base = 2, na.fill = 0) {
  model <- match.arg(model)
  
  if (!is.list(seqdef_res) || is.null(seqdef_res$seqdef)) {
    stop("seqdef_res must be a list containing a '$seqdef' element.")
  }
  
  s_scores <- seqdef_res$seqdef
  target_taxa <- names(s_scores)
  
  if (is.data.frame(trait_values)) {
    t_vec <- setNames(trait_values[[2]], as.character(trait_values[[1]]))
  } else {
    t_vec <- trait_values
  }
  if (is.null(names(t_vec))) {
    stop("trait_values must be a named vector (names matching taxa in seqdef_res$seqdef) ",
         "or a data frame with taxon names in column 1.")
  }
  
  missing_taxa <- setdiff(target_taxa, names(t_vec))
  if (length(missing_taxa) > 0) {
    warning(sprintf("%d taxa in seqdef_res have no trait value (e.g. %s); filled with na.fill = %s.",
                    length(missing_taxa),
                    paste(missing_taxa[seq_len(min(5, length(missing_taxa)))], collapse = ", "),
                    format(na.fill)))
  }
  
  aligned_traits <- t_vec[target_taxa]
  names(aligned_traits) <- target_taxa
  aligned_traits[is.na(aligned_traits)] <- na.fill
  
  if (model == "exponential") {
    priority <- s_scores * (base ^ aligned_traits)
  } else if (model == "linear") {
    priority <- s_scores * aligned_traits
  } else if (model == "additive") {
    priority <- s_scores + aligned_traits
  }
  
  return(priority)
}