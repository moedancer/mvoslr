#' Constructor for class "mvoslr_result" which is used to save and pass results
#'
#' @param raw_process Values of multivariate process at analysis dates
#' @param covariation_matrices Estimated covariance matrix for multivariate process at analysis date
#' @param multivariate_test_statistics Stagewise multivariate, standardised test-statistics
#' @param univariate_test_statistic Stagewise univariate test-statistics
#' @param stagewise_p_values Stagewise p-values
#' @param rejection_stage Number of stage at which the null hypothesis was rejected (is "NA" if null hypothesis could not be rejected)
#'
#' @return Object of class "mvoslr_result"
#'
#' @keywords internal
new_mvoslr_enhanced_result <- function(raw_process = array(),
                              stagewise_p_values = matrix(),
                              rejection_stage = integer(),
                              remaining_analyses = integer(),
                              vector_norm = character(),
                              variable_parameter = character(),
                              parameter_values = numeric(),
                              follow_up_fixed = logical()){

  stopifnot(is.array(raw_process))
  stopifnot(is.matrix(stagewise_p_values))
  stopifnot(is.integer(rejection_stage))
  stopifnot(is.integer(remaining_analyses))
  stopifnot(is.character(vector_norm))
  stopifnot(is.character(variable_parameter))
  stopifnot(is.numeric(parameter_values))
  stopifnot(is.logical(follow_up_fixed))

  dimname_template <- paste(variable_parameter, "=", parameter_values)

  dimnames(raw_process)[[3]] <- colnames(stagewise_p_values) <- names(rejection_stage) <- dimname_template

  result <- structure(list(raw_process = raw_process,
                           stagewise_p_values = stagewise_p_values,
                           rejection_stage = rejection_stage),
                      class = "mvoslr_enhanced_result",
                      remaining_analyses = remaining_analyses,
                      vector_norm = vector_norm,
                      variable_parameter = variable_parameter,
                      parameter_values = parameter_values,
                      follow_up_fixed = follow_up_fixed)

  return(result)

}

#' Print method for an object of class "mvoslr_enhanced_result"
#'
#' @param result Object of class "mvoslr_enhanced_result"
#'
#' @return Returns the object itself
#'
#' @keywords internal
print.mvoslr_enhanced_result <- function(result){

  # Only show (interim) test decisions here
  cat("Results of multivariate one-sample log-rank test with varying ", attributes(result)$variable_parameter, ".\n", sep = "")
  cat("Parameter assumes the following values:\n")
  print(attributes(result)$parameter_values)

  if(attributes(result)$variable_parameter == "a" & attributes(result)$fixed_follow_up){
    cat("\nDuration of follow-up period was held fixed.")
  } else if(attributes(result)$variable_parameter == "a" & attributes(result)$fixed_follow_up){
    cat("\nOverall duration of trial was held fixed,")
  }

  cat("\nFor these different parameters, the null hypothesis could be rejected in the follwoing stages (NA = acceptance of null hypothesis):\n")
  print(result$rejection_stage)

  return(result)

}

#' Summary method for an object of class "mvoslr_result"
#'
#' @param result Object of class "mvoslr_result"
#'
#' @return Returns the object itself
#'
#' @keywords internal
summary.mvoslr_enhanced_result <- function(result){

  # Show (interim) test decisions here
  cat("Results of multivariate one-sample log-rank test with varying ", attributes(result)$variable_parameter, ".\n", sep = "")
  cat("Parameter assumes the following values:\n")
  print(attributes(result)$parameter_values)

  if(attributes(result)$variable_parameter == "a" & attributes(result)$fixed_follow_up){
    cat("\nDuration of follow-up period was held fixed.")
  } else if(attributes(result)$variable_parameter == "a" & attributes(result)$fixed_follow_up){
    cat("\nOverall duration of trial was held fixed,")
  }

  cat("\nFor these different parameters, the null hypothesis could be rejected in the follwoing stages (NA = acceptance of null hypothesis):\n")
  print(result$rejection_stage)

  # Show values of multivariate process at analysis dates so far
  cat("Values of multivariate process at analysis dates for each parameter:\n")
  print(result$raw_process)

  # Show stagewise p-values at analysis dates so far
  cat("Stagewise p-values: (computed with", attributes(result)$vector_norm, "norm)\n")
  print(result$stagewise_p_values)

  return(result)

}
