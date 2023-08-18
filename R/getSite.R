#' Extract Phosphorylation Sites from Modified Peptide Sequence
#'
#' This function extracts the positions of phosphorylation sites from a modified peptide sequence.
#' It handles various modifications and returns the positions of the phosphorylation sites.
#'
#' Clean Input Sequence: The function removes specific modifications from the input sequence, such as "Carbamidomethyl (C)", "Acetyl (Protein N-term)", and "Oxidation (M)".
#' Split by Phosphorylation Sites: The cleaned sequence is split at each "[Phospho (STY)]" site, and the underscores are removed.
#' Calculate Positions: The function calculates the cumulative sum of the character lengths of the split strings, which represents the positions of the phosphorylation sites.
#' Return Positions: If only one phosphorylation site is found, the function returns NA. Otherwise, it returns the positions as a semicolon-separated string.
#' 
#' @param modseq A character string representing the modified peptide sequence.
#' @return A character string representing the positions of the phosphorylation sites, or NA if only one site is found.
#' @examples
#' getSite("_HSPAPPPDPGFPAPS[Phospho (STY)]PPPADSPSEGFS[Phospho (STY)]LK_")
#' getSite("_SKS[Phospho (STY)]PPKS[Phospho (STY)]PEEEGAVSS[Phospho (STY)]_")
#' @export
getSite <- function(modseq){
  #example: getSite("_HSPAPPPDPGFPAPS[Phospho (STY)]PPPADSPSEGFS[Phospho (STY)]LK_"
  #example" getSite("_SKS[Phospho (STY)]PPKS[Phospho (STY)]PEEEGAVSS[Phospho (STY)]_")
  
  cleaned=gsub("\\[Carbamidomethyl \\(C\\)\\]", '',
               gsub("\\[Acetyl \\(Protein N-term\\)\\]", '',
                    gsub("\\[Oxidation \\(M\\)\\]", '', modseq))) 
  
  #gsub('_', '', modseq)))
  
  x=cleaned %>%
    strsplit(., '\\[Phospho \\(STY\\)\\]') %>% 
    unlist %>%
    gsub('_', '', .) %>%
    nchar() %>%
    cumsum() 
  
  if(length(x)==1){
    return(NA)
  } else{
    return(paste0( x[1:length(x)-1], collapse=';'))
  }
}