#' Normalize Phosphorylation Sites Data
#'
#' This function performs various preprocessing and normalization steps on a dataset containing
#' information about phosphorylation sites. It reads the data, applies transformations, filters,
#' imputes missing values, and scales the data.
#'
#' @param sites A data frame containing phosphorylation sites information.
#' @param value_cutoff A numeric value representing the cutoff for log2 transformation (default is 100).
#' @param path An optional character string representing the path to read the data from.
#' @return A PhosphoExperiment object containing the normalized and processed data.
#' @export
normalizeData <- function(sites, value_cutoff = 100, path) {
  # Read data if path is provided
  if (!missing(path)) {
    sites <- data.table::fread(path)
    colnames(d)[1] <- 'id'
  }
  
  # Pivot and transform values
  d <- sites %>%
    tidyr::pivot_longer(!id, names_to = 'sample') %>%
    dplyr::mutate(value = ifelse(is.infinite(log2(value)) | value < value_cutoff, NA, log2(value))) %>%
    tidyr::pivot_wider(id_cols = id, names_from = sample, values_from = value)
  
  # Convert to matrix
  mat <- d %>%
    .[, -1] %>%
    as.matrix()
  rownames(mat) <- d %>%
    data.frame() %>%
    .[, 1]
  
  # Create PhosphoExperiment object
  ppe <- PhosphoExperiment(assays = list(Quantification = mat))
  
  # Filter, impute, and scale data
  grps <- stringr::str_sub(string = colnames(ppe), end = -4)
  ppe_filtered <- selectGrps(ppe, grps, percent = 0.5, n = 1)
  set.seed(123)
  ppe_imputed <- scImpute(ppe_filtered, 0.5, grps)[, colnames(ppe_filtered)]
  set.seed(123)
  ppe_imputed <- tImpute(ppe_imputed, assay = "imputed")
  ppe_imputed <- medianScaling(ppe_imputed, grps = F, scale = F, assay = 'imputed')
  
  # Return the result
  return(ppe_imputed)
}