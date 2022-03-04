#' Constructor for class "reference_model"
#'
#' @param transition_matrix Matrix of transitions between states as in mstate package
#' @param intensities Cumulative hazard functions for transitions in this model
#' @param type Reference multi-state model is either Markov (\code{model_type = "M"}) or Semi-Markov (\code{model_type = "SM"})
#' @param parametric Cumulative intensity functions can be parametric with Weibull form
#' @param parameters If cumulative intensity functions are parametric, parameters can be supplied here with one column for each transition, shape parameters in first and scale parameters in second row
#'
#' @return Object of class "reference model" to be used in planning and execution of multivariate one-sample log-rank test
#'
#' @export
#'
#' @examples
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' model_type_example <- "SM"
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cumhaz_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' parametric_example <- TRUE
#' parameters_example <- matrix(c(1.1, 1, 1.2, 1, 0.9, 1), nrow = 2,
#'                              dimnames = list(c("shape", "scale"), 1:3))
#' model_example <- new_reference_model(transition_matrix = tmat_example,
#'                                      intensities = cumhaz_example,
#'                                      type = model_type_example, parametric = parametric_example,
#'                                      parameters = parameters_example)
new_reference_model <- function(transition_matrix = matrix(), intensities = list(), type = character(),
                                parametric = FALSE, parameters = matrix()){

  stopifnot(is.matrix(transition_matrix))
  stopifnot(is.list(intensities))
  stopifnot(is.character(type))
  stopifnot(is.logical(parametric))
  stopifnot(is.matrix(parameters))

  model <- structure(list(transition_matrix = transition_matrix,
                          intensities = intensities,
                          parameters = parameters),
                     class = "reference_model",
                     type = type,
                     parametric = parametric)

  return(model)

}

#' Validator of the class "reference_model"
#'
#' @param model Object of class "reference_model"
#'
#' @return Object itself is returned if all checks are passed
#'
#' @export
#'
#' @examples
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' model_type_example <- "SM"
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cumhaz_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' parametric_example <- TRUE
#' parameters_example <- matrix(c(1.1, 1, 1.2, 1, 0.9, 1), nrow = 2,
#'                              dimnames = list(c("shape", "scale"), 1:3))
#' model_example <- new_reference_model(transition_matrix = tmat_example,
#'                                      intensities = cumhaz_example,
#'                                      type = model_type_example, parametric = parametric_example,
#'                                      parameters = parameters_example)
#' validate_reference_model(model_example)
validate_reference_model <- function(model){

  transmat <- model$transition_matrix
  matrix_size <- dim(transmat)

  if(matrix_size[1] != matrix_size[2]){
    stop("Transition matrix needs to be quadratic",
         call. = FALSE)
  }

  if(!all(model$transition_matrix) %in% c(NA, 1:dim(model$transition_matrix)[1]^2)){
    stop("Values of transition matrix need to be between 1 and the size of the matrix or NA",
         call. = FALSE)
  }

  if(!all(!duplicated(transmat[!is.na(transmat)]))){
    stop("Each value in transition matrix which is not NA has to be unique",
         call. = FALSE)
  }

  num_trans <- length(unique(as.vector(transmat[which(!is.na(transmat))])))
  intensities <- model$intensities

  if(length(intensities) != num_trans){
    stop("An intensity function for each transition is required",
         call. = FALSE)
  }

  if(! (attributes(model)$type %in% c("SM", "M")) ){
    stop("Model type needs to be either Markov (\"M\") or Semi-Markov (\"SM\")",
         .call = FALSE)
  }

  if(!all(unlist(lapply(intensities, is.function)))){
    stop("Elements of list of intensities need to be functions",
         .call = FALSE)
  }

  if(!all(unlist(lapply(intensities, function(fun) length(formals(fun)))) == 1)){
    stop("Each intensity has to be a univariate function",
         .call = FALSE)
  }

  if(attributes(model)$parametric){

    if(!identical(dim(model$parameters), as.integer(c(2, num_trans)))){
      stop("For a parametric model, a shape and a scale parameter for each transition is required",
           .call = FALSE)
    }

    if(!all(model$parameters > 0)){
      stop("Only positive parameter values are allowed",
           .call = FALSE)
    }

  }

  model

}

#' Helper function to create object of class "reference_model" from parameters of transition intensity functions with Weibull form.
#'
#' @param transition_matrix Matrix of transitions between states as in mstate package
#' @param type Reference multi-state model is either Markov (\code{model_type = "M"}) or Semi-Markov (\code{model_type = "SM"})
#' @param parameters Matrix of parameters for Weibull-shaped transition intensities with one column for each transition, shape parameters in first and scale parameters in second row
#'
#' @return Object of class "reference model" to be used in planning and execution of multivariate one-sample log-rank test
#' @export
#'
#' @examples
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' model_type_example <- "SM"
#' parameters_example <- matrix(c(1.1, 1, 1.2, 1, 0.9, 1), nrow = 2,
#'                              dimnames = list(c("shape", "scale"), 1:3))
#' model_example <- reference_model_weibull(transition_matrix = tmat_example,
#'                                          type = model_type_example,
#'                                          parameters = parameters_example)
reference_model_weibull <- function(transition_matrix = matrix(), type = character(),
                                    parameters = matrix()){

  stopifnot(is.matrix(transition_matrix))
  stopifnot(is.character(type))
  stopifnot(is.matrix(parameters))

  transmat <- transition_matrix
  matrix_size <- dim(transmat)

  parameters_list <- lapply(1:ncol(parameters),
                                     function(col_number) parameters[, col_number])
  intensities <- lapply(parameters_list,
                        FUN = get_cum_haz_fct_weibull)

  model <- structure(list(transition_matrix = transition_matrix,
                          intensities = intensities,
                          parameters = parameters),
                     class = "reference_model",
                     type = type,
                     parametric = TRUE)

  if(class(validate_reference_model(model)) == "reference_model"){
    return(model)
  } else {
    stop("Reference model could not be created due to failed validation",
         .call = FALSE)
  }

}


#' Print method for an object of class "reference_model"
#'
#' @param model Object of class "reference_model"
#'
#' @return Returns the object itself
#'
#' @keywords internal
print.reference_model <- function(model){

  long_form_parametric <- ifelse(attributes(model)$parametric, "Parametric", "Non-parametric")
  long_form_type <- ifelse(attributes(model)$type == "SM", "Semi-Markov", "Markov")

  num_states <- dim(model$transition_matrix)[1]
  num_trans <- length(model$intensities)

  cat(long_form_parametric, long_form_type, "model with", num_states, "states and", num_trans, "transitions" )

  invisible(model)

}

#' Summary method for an object of class "reference_model"
#'
#' @param model Object of class "reference_model"
#'
#' @return Returns the object itself
#'
#' @keywords internal
summary.reference_model <- function(model){

  long_form_type <- ifelse(attributes(model)$type == "SM", "Semi-Markov", "Markov")

  cat(long_form_type, "model with transitions\n")

  print(mstate::to.trans2(model$transition_matrix))

  if(!attributes(model)$parametric){
    cat("Model is not specified as parametric. Cumulative intensities can be looked at with the \"plot\" function.\n")
  } else {
    cat("Transition intensities have Weibull shape with the following parameters.\nCumulative intensities can be looked at with the \"plot\" function.\n")
    print(model$parameters)
  }

  invisible(model)

}

#' Plot method for an object of class "reference_model". Plots cumulative intensity functions of all transitions.
#'
#' @param model Object of class "reference_model"
#'
#' @return Returns the object itself
#'
#' @import graphics
#'
#' @keywords internal
plot.reference_model <- function(model, ...){

  num_trans <- length(model$intensities)

  plot(model$intensities[[1]], main = "Cumulative intensity functions of all transitions",
       xlab = "Time", ylab = "Cumulative intensity", ...)

  col_count <- 1

  if(length(model$intensities) > 1){

    for(f_index in 2:length(model$intensities)){
      f <- model$intensities[[f_index]]
      curve(f, add = TRUE, col = f_index, ...)
    }

  }

  legend("topleft", legend = 1:num_trans, col = 1:num_trans, lty = 1)

  invisible(model)

}
