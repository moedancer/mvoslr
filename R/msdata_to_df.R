#' Transform data of type msdata from mstate package to a data frame
#'
#' @param msdata
#'
#' @return data frame with column names corresponding to the columns names of msdata objects
msdata_to_df <- function(msdata){

  ms_names <- names(msdata)
  attributes(msdata) <- NULL

  ms_df <- data.frame(msdata)
  colnames(ms_df) <- ms_names

  return(ms_df)

}
