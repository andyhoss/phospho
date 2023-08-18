#' Retrieve Top Hits from Fitted Linear Model
#'
#' This function retrieves the top hits from a fitted linear model object based on specified
#' contrasts, p-value threshold, and log fold change (LFC) threshold. It is commonly used in
#' differential expression analysis.
#'
#' @param fit An object representing the fitted linear model (e.g., from lmFit).
#' @param contrasts A character string or expression representing the contrasts to be used in the analysis.
#' @param lfc A numeric value representing the log fold change threshold for filtering the results.
#' @return A data frame containing the top hits based on the specified criteria.
#' @export
getTopHits <- function(fit, contrasts, lfc) {
  # Retrieve top hits using specified contrasts, p-value threshold, and LFC threshold
  tophits <- topTable(fit, coef = contrasts, p.value = 0.05, lfc = lfc, n = Inf)
  
  # Return the result
  return(tophits)
}
