#' Perform Linear Model Fitting and Hypothesis Testing
#'
#' This function performs a linear model fitting and hypothesis testing on a given data object
#' using specified contrasts. It is designed for analyzing differential expression in proteomics data.
#'
#' @param ppe A specialized data object containing proteomics information.
#' @param contrasts A character string or expression representing the contrasts to be used in the analysis.
#' @return An object representing the fitted model with contrasts and Empirical Bayes smoothing applied.
#' @export
testPPE <- function(ppe, contrasts){
  grps <- stringr::str_sub(string = colnames(ppe), end=-4)
  design <- model.matrix(~ grps - 1)
  fit <- lmFit(ppe@assays@data$scaled, design)
  contrast.matrix <- makeContrasts(contrasts=contrasts, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  return(fit2)
}