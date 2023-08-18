#' Extract Window of Peptide Sequence
#'
#' This function extracts a specific window of a peptide sequence based on a given position
#' and window size. It handles edge cases by padding with underscores.
#'
#'Calculate Window Bounds: The function calculates the left and right boundaries of the window based on the provided position and window size.
#' Handle Missing Values: If either the position or the peptide is missing (NA), the function sets the region to NA.
#' Extract Window: If the left and right sides of the window are within the peptide bounds, the function extracts the corresponding substring.
#' Handle Left Boundary: If the left side of the window is less than 1, the function pads the region with underscores.
#' Handle Right Boundary: If the right side of the window is greater than the length of the peptide, the function pads the region with underscores.
#' Return Region: The function returns the extracted region.
#' 
#' @param pos An integer representing the central position of the window.
#' @param peptide A character string representing the peptide sequence.
#' @param window_size An integer representing the size of the window (default is 15).
#' @return A character string representing the extracted window of the peptide sequence.
#' @export
getWindow <- function(pos, peptide, window_size=15){
  left_side=pos-((window_size-1)/2)
  right_side=pos+((window_size-1)/2)
  
  if(is.na(pos) | is.na(peptide)){
    region=NA
  }
  if(left_side >=1 & right_side <= nchar(peptide)){
    region=substring(first=left_side, last=right_side, peptide)
  }
  
  if(left_side < 1){ 
    region = paste0(c(rep('_', abs(left_side-1)), substring(first=left_side, last=right_side, peptide)), collapse='')
  }
  
  if(right_side > nchar(peptide)){
    temp = substring(text=peptide, first=pos-7)
    region = paste0(c(temp, rep(x='_', (window_size-nchar(temp)))), collapse='')
  }
  return(region)
}