#' Helper function to create cumulative hazard function for exponential distribution
#'
#' @param rate Rate of exponential distribution
#'
#' @return conditional hazard function
#'
#' @keywords internal
#'
#' @examples
#' rate <- 1
#' reference_ch_fct <- get_cum_haz_fct_exp(rate)
get_cum_haz_fct_exp <- function(rate){

  rate <- unname(rate)
  ch_fct <- function(t) rate*t
  return(ch_fct)

}

#' Helper function to create cumulative hazard function for exponential distribution
#'
#' @param parameters Vector of length two containing shape and scale parameter
#'
#' @return Conditional hazard function
#'
#' @keywords internal
#'
#' @examples
#' shape <- 1.5
#' scale <- 1
#' reference_ch_fct <- get_cum_haz_fct_weibull(c(shape, scale))
get_cum_haz_fct_weibull <- function(parameters){

  parameters <- unname(parameters)
  ch_fct <- function(t) (t/parameters[2])^parameters[1]
  return(ch_fct)

}
