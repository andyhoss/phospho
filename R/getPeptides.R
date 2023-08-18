#' Retrieve Peptide Sequences
#'
#' This function retrieves peptide sequences based on provided accessions. It uses various sources
#' including BioMart and custom functions to get the sequences, and handles missing values.
#' 
#' Read Accessions: The function reads a file containing protein accessions, splitting and unnesting them as needed.
#' Query Peptide Sequences: It queries the BioMart database to get peptide sequences for the given accessions using the query_biomart function and a specific filter attribute.
#' Load Peptide Sequences: The function loads previously saved peptide sequences from an RData file (this part is commented out in the code).
#' Find Missing Peptides: It finds the missing peptide sequences and attempts to recover them.
#' Query Other Accessions: The function performs additional queries using different filter attributes to get more peptide sequences.
#' Process Peptide Sequences: It processes the retrieved peptide sequences, removing any trailing stars.
#' Retrieve Remaining Missing Sequences: The function uses a loop to try to retrieve any remaining missing sequences using the get_protein_sequence function.
#' Combine Final Sequences: It combines all the retrieved sequences into a final result.
#' Return Final Sequences: The function returns the final peptide sequences.
#' 
#' @param path_to_raw_data A character string representing the path to the raw data file containing protein accessions.
#' @return A data frame containing the final peptide sequences.
#' @export
getPeptides <- function(path_to_raw_data){
  accs <- data.table::fread(path_to_raw_data) %>% 
    dplyr::select(PG.ProteinAccessions) %>%
    dplyr::mutate(accessions = strsplit(PG.ProteinAccessions, ';')) %>% 
    tidyr::unnest(accessions) %>% 
    dplyr::distinct()
  
  pep_seq <- phospho::query_biomart(data=accs$accessions, filter_attribute = "uniprotswissprot")
  
  #save(pep_seq, file='pep_seq.RData')
  load('pep_seq.RData')
  #"uniprot_gn_symbol", 
  #"uniprotsptrembl",
  
  
  missing_pgs <- accs %>% 
    dplyr::left_join(pep_seq, by=c('accessions'='uniprotswissprot')) %>% 
    dplyr::filter(is.na(peptide)) %>% 
    pull(PG.ProteinAccessions) %>% 
    unique()
  
  recovered_pg <- accs %>% 
    dplyr::left_join(pep_seq, by=c('accessions'='uniprotswissprot')) %>% 
    dplyr::filter(PG.ProteinAccessions %in% missing_pgs) %>% 
    dplyr::filter(!is.na(peptide)) %>% 
    pull(PG.ProteinAccessions) %>%
    unique()
  
  missing_acc <- accs %>% 
    dplyr::left_join(pep_seq, by=c('accessions'='uniprotswissprot')) %>% 
    dplyr::filter(PG.ProteinAccessions %in% setdiff(missing_pgs, recovered_pg)) %>% 
    dplyr::pull(accessions) %>% 
    unique()
  
  other_acc1 <- phospho::query_biomart(data=missing_acc, filter_attribute = "uniprotsptrembl")
  other_acc2 <- phospho::query_biomart(data=missing_acc, filter_attribute = "uniprotswissprot")
  
  temp_acc <- dplyr::bind_rows(pep_seq %>% 
                                 dplyr::rename(accessions = uniprotswissprot), 
                               other_acc1 %>% 
                                 dplyr::rename(accessions = uniprotsptrembl),
                               other_acc2 %>%
                                 dplyr::rename(accessions = uniprotswissprot)) %>%
    dplyr::select(accessions, peptide) %>%
    dplyr::mutate(peptide = stringr::str_sub(peptide, end=-2)) #remove star from end of peptide seq
  
  still_missing <- setdiff(missing_acc, temp_acc$accessions)
  
  other_acc3 <- list()
  for(i in 1:length(still_missing)){
    other_acc3[[i]] <- data.frame(accessions = still_missing[i],
                                  peptide = tryCatch( expr = { 
                                    phospho::get_protein_sequence(still_missing[i]) }, 
                                    error = function(e) { NA })
    )
  }
  
  final_seq <- temp_acc %>%
    dplyr::bind_rows(
      dplyr::bind_rows(other_acc3) %>% 
        dplyr::filter(!is.na(peptide))
    )
  return(final_seq)
}