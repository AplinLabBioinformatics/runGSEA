#' Create .cls File
#'
#' Creates a .cls file for running GSEA
#' 
#' @param annotDF A dataframe object with corresponding sample information and a column with sample groups
#' @param groupingColName The column with group names for the samples
#' @param filenameOut Destination of file
#' @details
#' The sample annotation dataframe with a column - default to first - that indicates the group a sample belongs. 
#' 
#' @return silently returns NULL
#' @export

# Adapted from:
# https://github.com/BenaroyaResearch/geneSetTools/blob/master/R/write_cls.R

write_cls <- function(annotDF, filenameOut, groupingColName = 1){
  nSamples = dim(annotDF)[1]
  nClasses = length(unique(annotDF[,groupingColName]))
  cat(sprintf("%d %d 1\n", nSamples, nClasses), file = filenameOut)
  ucl_1 = paste(unique(annotDF[,groupingColName]), collapse = " ")
  cat(paste0("# ", ucl_1), sep = "\n", append = TRUE, file = filenameOut)
  
  ucl_A = paste(annotDF[,groupingColName], collapse = " ")
  cat(ucl_A, sep = "\n", append = TRUE, file = filenameOut)
}

