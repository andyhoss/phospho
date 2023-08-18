#' Retrieve Significant Motifs from Fitted Linear Model
#'
#' This function retrieves significant motifs from a fitted linear model object based on specified
#' contrasts, log fold change (LFC) threshold, and sites of interest. It performs a motif search
#' using the specified foreground and background sequences and returns the identified motifs.
#'
#' @param fit An object representing the fitted linear model (e.g., from lmFit).
#' @param sites A character string representing the sites of interest (default is "STY").
#' @param data A data frame containing information about the sequences, including 'id' and 'window' columns.
#' @param contrasts A character string or expression representing the contrasts to be used in the analysis.
#' @param LFC A numeric value representing the log fold change threshold for filtering the results (default is 1.5).
#' @return An object representing the identified motifs based on the specified criteria.
#' @export
getMotif <- function(fit, sites = "STY", data, contrasts, LFC = 1.5) {
  # Retrieve top hits
  top <- getTopHits(fit, contrasts = contrasts, lfc = LFC)
  
  # Get paths to sample files
  bg.path = system.file("extdata", "bg-data-serine.txt", package = "rmotifx")
  
  # Read in background, serine centered sequences
  bg.seqs = readLines(bg.path)
  
  # Prepare foreground sequences
  fg.seqs = data %>%
    dplyr::select(id, window) %>%
    dplyr::distinct() %>%
    dplyr::filter(id %in% rownames(top %>% dplyr::filter(logFC < 0))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(consensus = ifelse(nchar(window) == 15, window,
                                     paste0(seqinr::consensus(stringr::str_split(unlist(stringr::str_split(window, ',')), '', simplify = T)), collapse = ''))) %>%
    dplyr::pull(consensus) %>%
    unique()
  
  # Run motif search
  mot = rmotifx::motifx(fg.seqs, bg.seqs, central.res = sites, min.seqs = 20, pval.cutoff = 1e-6)
  
  # Return the result
  return(mot)
}
