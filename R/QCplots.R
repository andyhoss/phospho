#' Create Quality Control (QC) Plots
#'
#' This function creates quality control (QC) plots for a given data object. It generates
#' both dendrogram and principal component analysis (PCA) plots for the quantification and
#' scaled assays within the data object.
#'
#' @param ppe A specialized data object containing proteomics information.
#' @return A ggplot object representing the arranged QC plots (dendrogram and PCA).
#' @export
QCplots <- function(ppe) {
  # Extract group information from column names
  grps <- stringr::str_sub(string = colnames(ppe), end = -4)
  
  # Create dendrogram plots for quantification and scaled assays
  p1 <- plotQC(SummarizedExperiment::assay(ppe, "Quantification"),
               grps = grps,
               labels = colnames(ppe),
               panel = "dendrogram")
  
  p2 <- plotQC(SummarizedExperiment::assay(ppe, "scaled"),
               grps = grps,
               labels = colnames(ppe),
               panel = "dendrogram")
  
  # Create PCA plots for quantification and scaled assays
  p3 <- plotQC(SummarizedExperiment::assay(ppe, "Quantification"),
               grps = grps,
               labels = colnames(ppe),
               panel = "pca")
  
  p4 <- plotQC(SummarizedExperiment::assay(ppe, "scaled"),
               grps = grps,
               labels = colnames(ppe),
               panel = "pca")
  
  # Arrange and display plots
  ggpubr::ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
}
