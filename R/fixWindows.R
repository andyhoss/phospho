#' Fix Missing Windows in Peptide Fragments
#'
#' This function processes a data frame containing peptide fragments and attempts to fix missing window
#' information. It identifies missing windows and iteratively calls other functions to fix them.
#'
#' Identify Missing Windows: The function identifies peptide fragments with missing window information and joins them with the original fragments data frame.
#' Fix Missing Windows (First Pass): The function iteratively calls the fixResidue function on each missing window to attempt to fix it. The fixed windows are then combined into a data frame.
#' Handle Still Missing Windows: The function identifies windows that are still missing after the first pass and attempts to fix them by retrieving protein sequences using the get_protein_sequence function and extracting windows using the getWindow function.
#' Determine Matches: The function determines matches by comparing various substrings of the stripped peptide sequence and window.
#' Fix Missing Windows (Second Pass): A second pass of fixing is performed, similar to the first pass, on the windows that are still missing.
#' Combine Fixed Windows: The function combines the fixed windows from both passes and filters those that match.
#' Combine with Original Fragments: The function combines the fixed windows with the original fragments that already had matching windows.
#' Return Result: The final result is returned, containing the fixed windows along with other associated information.
#'
#' @param fragments A data frame containing peptide fragments with associated information.
#' @return A data frame containing the fixed windows along with other associated information.
#' @export
fixWindows <- function(fragments){
  missing_windows <- fragments %>% 
    dplyr::group_by(PG.ProteinAccessions, EG.ProteinPTMLocations) %>% 
    dplyr::summarize(match = max(match)) %>% 
    dplyr::filter(match==0) %>%
    dplyr::left_join(fragments %>% 
                       dplyr::select(-match), by=c('PG.ProteinAccessions', 'EG.ProteinPTMLocations'))
  
  
  fixed_windows <- list()
  for(i in 1:nrow(missing_windows)){
    fixed_windows[[i]] <- phospho::fixResidue(fragment=missing_windows[i,] %>% pull(PEP.StrippedSequence),
                                     peptide =missing_windows[i,] %>% pull(peptide),
                                     site = missing_windows[i,] %>% pull(sites))
    print(i)
  }
  fixed_windows <- dplyr::bind_rows(fixed_windows)
  
  still_missing <- fragments %>% 
    dplyr::inner_join(fixed_windows %>% 
                        dplyr::group_by(PEP.StrippedSequence, sites) %>% 
                        dplyr::summarise(match = max(match)) %>% 
                        dplyr::filter(match == 0)  %>% 
                        dplyr::select(-match)) 
  
  swissprot_selected <- FALSE
  new_pep <- still_missing %>%
    dplyr::select(accessions) %>%
    dplyr::distinct() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(peptide = phospho::get_protein_sequence(accessions))
  
  still_missing <- still_missing %>% 
    dplyr::select(-c(peptide, window, match)) %>% 
    dplyr::left_join(new_pep, by='accessions') %>%
    dplyr::rowwise() %>%
    dplyr::mutate(window = phospho::getWindow(pos=pos, peptide=peptide))%>% 
    dplyr::mutate(match = ifelse(nchar(substring(PEP.StrippedSequence, as.numeric(sites)-2, as.numeric(sites)+2))==5 &
                                   substring(PEP.StrippedSequence, as.numeric(sites)-2, as.numeric(sites)+2) == substring(window,6,10), 1, 
                                 ifelse(substring(PEP.StrippedSequence, as.numeric(sites)-1, as.numeric(sites)+1) == substring(window,7,9), 1, 
                                        ifelse(sites ==1 & substring(PEP.StrippedSequence, as.numeric(sites), as.numeric(sites)+2) == substring(window,8,10), 1, 
                                               ifelse(pos == nchar(peptide) & substring(PEP.StrippedSequence, as.numeric(sites)-2, as.numeric(sites)) == substring(window,6,8), 1,
                                                      0))))) %>%
    dplyr::ungroup()
  
  missing_windows2 <- still_missing %>% 
    dplyr::group_by(PG.ProteinAccessions, EG.ProteinPTMLocations) %>% 
    dplyr::summarize(match = max(match)) %>% 
    dplyr::filter(match==0) %>%
    dplyr::left_join(still_missing %>% 
                       dplyr::select(-match), by=c('PG.ProteinAccessions', 'EG.ProteinPTMLocations'))
  
  fixed_windows2 <- list()
  for(i in 1:nrow(missing_windows2)){
    fixed_windows2[[i]] <- phospho::fixResidue(fragment=missing_windows2[i,] %>% pull(PEP.StrippedSequence),
                                      peptide =missing_windows2[i,] %>% pull(peptide),
                                      site = missing_windows2[i,] %>% pull(sites))
    print(i)
  }
  fixed_windows2 <- dplyr::bind_rows(fixed_windows2) 
  
  final_fixed <- bind_rows(fixed_windows, fixed_windows2) %>% 
    dplyr::filter(match==1) %>% 
    select(-match) %>% 
    left_join(fragments %>% 
                dplyr::select(PG.ProteinAccessions, PEP.StrippedSequence, FG.LabeledSequence, EG.ProteinPTMLocations, accessions, sites, site_count, m), 
              by=c('PEP.StrippedSequence', 'sites')) 
  
  out <- final_fixed %>%
    dplyr::bind_rows(fragments %>% dplyr::filter(match==1) %>% dplyr::select(-match)) %>%
    dplyr::arrange(PG.ProteinAccessions, EG.ProteinPTMLocations) %>%
    dplyr::select(PG.ProteinAccessions, EG.ProteinPTMLocations, FG.LabeledSequence, PEP.StrippedSequence, 
                  sites, accessions, ptm, site_count, residue, m, pos, peptide, window) 
  
  return(out)
  
}