#' Create Distribution Plots
#'
#' This function generates a histogram plot for the distribution of values from a specified data object.
#' It provides options for using scaled values and aggregating results.
#'
#' @param ppe A specialized data object containing proteomics information.
#' @param aggregate A logical value indicating whether to aggregate the results (default is TRUE).
#' @param scaled A logical value indicating whether to use scaled values (default is TRUE).
#' @return A ggplot object representing the histogram plot.
#' @export
distPlots <- function(ppe, aggregate=TRUE, scaled=TRUE){
  if(scaled==TRUE){
    #distribution plot
    p1=ppe@assays@data@listData$scaled %>%
      data.frame() %>%
      tibble::rownames_to_column('site') %>% 
      tidyr::pivot_longer(!site, names_to ='sample') %>% 
      ggplot(., aes(x=value)) + 
      geom_histogram(bins=200) +
      #facet_grid(rows = vars(sample)) +
      theme(strip.text.y = element_text(angle = 0))
  } else{
    p1=ppe@assays@data@listData$Quantification %>%
      data.frame() %>%
      tibble::rownames_to_column('site') %>% 
      tidyr::pivot_longer(!site, names_to ='sample') %>% 
      ggplot(., aes(x=value)) + 
      geom_histogram(bins=200) +
      theme(strip.text.y = element_text(angle = 0))
  }
  if(aggregate==FALSE){
    p1 <- p1 +       
      facet_grid(rows = vars(sample)) 
  }
  return(p1)
}