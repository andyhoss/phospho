#' Retrieve Protein Sequence from Swiss-Prot
#'
#' This function retrieves a protein sequence from the Swiss-Prot database using a given accession number.
#' It ensures that the Swiss-Prot bank is selected before performing the query.
#'
#' @import seqinr
#' @param accession_number A character string representing the accession number of the protein.
#' @return A character string representing the protein sequence.
#' @export
get_protein_sequence <- function(accession_number) {
  global <- get("swissprot_selected", envir = globalenv())
  #swissprot_selected <- FALSE
  # Check if the Swiss-Prot bank is already chosen
  if (!global) {
    # Choose the Swiss-Prot server if not already chosen
    seqinr::choosebank("swissprot")
    assign("swissprot_selected", TRUE, envir = globalenv())
  }
  
  # Query the sequence using the accession number
  protein_query <- seqinr::query("list", paste("AC=", accession_number, sep=""))
  
  # Get the sequence
  protein_sequence <- phospho::getSequence(protein_query$req[[1]], as.string = TRUE)
  
  # Return the sequence
  return(protein_sequence)
}