#' Fix Residue within Peptide Fragment
#'
#' This function attempts to fix a residue within a peptide fragment by performing a local alignment
#' between the fragment and the full peptide sequence. It calculates the new position, residue, and window
#' and returns the fixed residue information.
#'
#' Perform Local Alignment: The function uses the pairwiseAlignment function from the Biostrings package to perform a local alignment between the peptide and the fragment.
#' Calculate New Position: The new position of the residue is calculated based on the alignment's start position and the provided site.
#' Extract New Residue: The new residue is extracted from the peptide sequence at the calculated position.
#' Extract New Window: The getWindow function is called to extract a window around the new position from the peptide sequence.
#' Determine Match: The function determines whether the fragment and the new window match based on various conditions.
#' Return Result: The function returns a data frame containing the fixed residue information, including the fragment, peptide, site, residue, position, window, and match status.
#' 
#' @param fragment A character string representing the peptide fragment.
#' @param peptide A character string representing the full peptide sequence.
#' @param site An integer representing the site position of the residue within the fragment.
#' @return A data frame containing the fixed residue information, including fragment, peptide, site, residue, position, window, and match status.
#' @export 
fixResidue <- function(fragment, peptide, site){
  x=Biostrings::pairwiseAlignment( peptide, fragment, type='local')
  
  #x@pattern %>% as.character() %>% substring(., site, site)
  
  new_pos = x@pattern@range@start + as.numeric(site) - 1 
  new_residue = paste0(substring(peptide, new_pos, new_pos), new_pos)
  new_window = getWindow(pos = new_pos, 
                         peptide = peptide, 
                         window_size = 15)
  
  match=ifelse(substring(fragment, as.numeric(site)-1, as.numeric(site)+1) == substring(new_window,7,9), 1, 
               ifelse(as.numeric(site) ==1 & substring(fragment, as.numeric(site), as.numeric(site)+2) == substring(new_window,8,10), 1, 
                      ifelse(new_pos == nchar(peptide) & substring(fragment, as.numeric(site)-2, as.numeric(site)) == substring(new_window,6,8), 1,
                             0)))
  
  return(data.frame(
    PEP.StrippedSequence=fragment,
    peptide = peptide,
    sites = site,
    residue = new_residue,
    pos = new_pos,
    window = new_window,
    match = match))
}