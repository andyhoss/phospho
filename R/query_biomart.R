#' Query BioMart by Dividing Data into Chunks
#'
#' This function queries the BioMart database by dividing the input data into chunks of 200 elements each.
#'
#' @param data A vector containing the data to be queried (e.g., UniProt accession IDs).
#' @param filter_attribute A filter attribute to be used in the BioMart query.
#' @return A data frame containing the results of the BioMart query.
#' \dontrun{
#'   data <- c('P51587', 'Q8N726') # Add more UniProt accession IDs as needed
#'   result <- query_biomart(data, filter_attribute = "uniprotswissprot")
#' }
#' @export
query_biomart <- function(data, filter_attribute) {
  ## Example usage
  ##data <- c('P51587', 'Q8N726') # Add more UniProt accession IDs as needed
  ##result <- query_biomart(data)
  
  # Determine the number of chunks
  num_chunks <- ceiling(length(data) / 200)
  print(paste("Start:", num_chunks, "batch to process"))
  
  # Connect to the desired mart (Ensembl)
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")
  
  # Use the desired dataset (for example, human genes)
  ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = mart)
  
  # Loop through each chunk and perform query
  results <- lapply(1:num_chunks, function(i) {
    # Select the data for this chunk
    chunk <- data[((i - 1) * 200 + 1):min(i * 200, length(data))]
    
    # Define the filters and values for the query
    attributes <- c("uniprotswissprot", "uniprotsptrembl","peptide")
    filters <- filter_attribute 
    values <- chunk
    
    # Query BioMart
    result <- biomaRt::getBM(attributes = attributes, filters = filters, values = values, mart = ensembl)
    print(paste("Processed: batch", i))
    return(result)
  })
  
  # Combine the results
  final_result <- do.call(rbind, results)
  print("Completed")
  return(final_result)
}