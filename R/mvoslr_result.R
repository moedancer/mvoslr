#' Constructor for class "mvoslr_result" which is used to save and display results of multivariate one-sample log-rank tests
#'
#' @param raw_process Values of multivariate process at analysis dates
#' @param covariance_matrices Estimated covariance matrix for multivariate process at analysis date
#' @param multivariate_test_statistics Stagewise multivariate, standardised test-statistics
#' @param univariate_test_statistic Stagewise univariate test-statistics
#' @param stagewise_p_values Stagewise p-values
#' @param rejection_stage Number of stage at which the null hypothesis was rejected (is \code{NA} if null hypothesis could not be rejected)
#'
#' @return Object of class "mvoslr_result"
#'
#' @keywords internal
new_mvoslr_result <- function(raw_process = matrix(),
                              covariance_matrices = array(),
                              multivariate_test_statistics = matrix(),
                              univariate_test_statistics = numeric(),
                              stagewise_p_values = numeric(),
                              rejection_stage = integer(),
                              remaining_analyses = integer(),
                              vector_norm = character()){

  stopifnot(is.matrix(raw_process))
  stopifnot(is.array(covariance_matrices))
  stopifnot(is.matrix(multivariate_test_statistics))
  stopifnot(is.numeric(univariate_test_statistics))
  stopifnot(is.numeric(stagewise_p_values))
  stopifnot(is.integer(rejection_stage) | is.na(rejection_stage))
  stopifnot(is.integer(remaining_analyses))
  stopifnot(is.character(vector_norm))

  # If events have no names in "raw_process", give them default names and apply those names to variance component
  num_events <- dim(raw_process)[1]
  default_event_names <- paste("Event", 1:num_events)
  if(is.null(dimnames(raw_process)[[1]])){
    dimnames(raw_process)[[1]] <- default_event_names
    dimnames(covariance_matrices)[[1]] <- dimnames(covariance_matrices)[[2]] <- default_event_names
    dimnames(multivariate_test_statistics)[[1]] <- default_event_names
  }

  # If analysis dates have no names in "raw_process", give them default names and apply those names to other components
  num_analyses <- dim(raw_process)[2]
  default_date_names <- paste("Analysis", 1:num_analyses)
  if(is.null(dimnames(raw_process)[[2]])){
    dimnames(raw_process)[[2]] <- default_date_names
    dimnames(covariance_matrices)[[3]] <- dimnames(multivariate_test_statistics)[[2]] <- default_date_names
    names(univariate_test_statistics) <- names(stagewise_p_values) <- default_date_names
  }

  result <- structure(list(raw_process = raw_process,
                           covariance_matrices = covariance_matrices,
                           multivariate_test_statistics = multivariate_test_statistics,
                           univariate_test_statistics = univariate_test_statistics,
                           stagewise_p_values = stagewise_p_values,
                           rejection_stage = rejection_stage),
                      class = "mvoslr_result",
                      remaining_analyses = remaining_analyses,
                      vector_norm = vector_norm)

  return(result)

}

#' Print method for an object of class "mvoslr_result"
#'
#' @param result Object of class "mvoslr_result"
#'
#' @return Returns the object itself
#'
#' @keywords internal
print.mvoslr_result <- function(result){

  # Only show (interim) test decision here
  if(is.na(result$rejection_stage) & result$remaining_analyses == 0){
    cat("The null hypothesis could not be rejected. H0 is accepted.")
  } else if(is.na(result$rejection_stage) & result$remaining_analyses > 0) {
    cat("The null hypothesis could not be rejected so far. Continue trial.")
  } else {
    cat("The null hypothesis can be rejected in stage ", result$rejection_stage, ".", sep = "")
  }

  invisible(result)

}

#' Summary method for an object of class "mvoslr_result"
#'
#' @param result Object of class "mvoslr_result"
#'
#' @return Returns the object itself
#'
#' @keywords internal
summary.mvoslr_result <- function(result){

  # Display (interim) test decision
  if(is.na(result$rejection_stage) & result$remaining_analyses == 0){
    cat("The null hypothesis could not be rejected. H0 is accepted.")
  } else if(is.na(result$rejection_stage) & result$remaining_analyses > 0) {
    cat("The null hypothesis could not be rejected so far. Continue trial.")
  } else {
    cat("The null hypothesis can be rejected in stage", result$rejection_stage)
  }

  # Show values of multivariate process at analysis dates so far
  cat("Values of multivariate process at analysis dates:\n")
  print(result$raw_process)

  # Show stagewise p-values at analysis dates so far
  cat("Stagewise p-values: (computed with", attributes(result)$vector_norm, "norm)\n")
  print(result$stagewise_p_values)

  return(result)

}
