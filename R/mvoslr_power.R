#' Constructor for class "mvoslr_power_object" which is used to save and pass results of power simulations.
#' If simulations for several values of a design parameter are conducted, results for each value are saved here.
#'
#' @param power Overall simulated power
#' @param rejection_stages Relative frequency of rejections of the null hypothesis in each stage
#' @param means Mean value of the multivariate process at analysis dates
#' @param variances Mean estimated covariance matrices of the multivariate process at analysis dates
#' @param runs Number of simulation runs
#' @param variable_parameter Varied parameter (if any)
#' @param follow_up_fixed Boolean variable indicating if the follow-up duration is held fixed when the accrual duration ("a") is varied
#' @param parameter_values Values of the parameter which is varied
#'
#' @return Object of class "mvoslr_result"
#'
#' @keywords internal
new_mvoslr_power_object <- function(power = numeric(),
                                    rejection_stages = matrix(),
                                    means = array(),
                                    variances = array(),
                                    runs = numeric(),
                                    variable_parameter = character(),
                                    follow_up_fixed = logical(),
                                    parameter_values = numeric()){

  stopifnot(is.numeric(power))
  stopifnot(is.matrix(rejection_stages))
  stopifnot(is.array(means))
  stopifnot(is.array(variances))
  stopifnot(is.numeric(runs))
  stopifnot(is.character(variable_parameter) | is.na(variable_parameter))
  stopifnot(is.numeric(parameter_values))
  stopifnot(is.logical(follow_up_fixed))

  # If there is a variable parameter, name results correspondingly
  if(length(variable_parameter) > 0 &
     !is.na(variable_parameter) &
     !is.null(variable_parameter)){

    dimname_template <- paste(variable_parameter, "=", parameter_values)
    names(power) <- colnames(rejection_stages) <- dimnames(means)[[3]] <- dimname_template

    # If variances are given, also name them correspondingly
    if(length(dim(variances)) == 4){
      dimnames(variance)[[4]] <- dimname_template
    }

  }

  # If events have no names in "means", give them default names and apply those names to variance component (if possible)
  num_events <- dim(means)[1]
  default_event_names <- paste("Event", 1:num_events)
  if(is.null(dimnames(means)[[1]])){
    dimnames(means)[[1]] <- default_event_names
    if(length(dim(variances)) %in% 3:4){
      dimnames(variances)[[1]] <- dimnames(variances)[[2]] <- default_event_names
    }
  }

  # If analysis dates have no names in "means", give them default names and apply those names to other components
  num_analyses <- dim(means)[2]
  default_date_names <- paste("Analysis", 1:num_analyses)
  if(is.null(dimnames(means)[[2]])){
    dimnames(means)[[2]] <- default_date_names
    if(length(dim(variances)) %in% 3:4){
      dimnames(variances)[[3]] <- default_event_names
    }
  }

  power_obj <- structure(list(power =  power,
                              rejection_stages = rejection_stages,
                              means = means,
                              variances = variances),
                         class = "mvoslr_power_object",
                         simulation_runs = runs,
                         variable_parameter = variable_parameter,
                         parameter_values = parameter_values,
                         follow_up_fixed = follow_up_fixed)

  return(power_obj)

}

#' Print method for an object of class "mvoslr_power_object"
#'
#' @param power_obj Object of class "mvoslr_power_object"
#'
#' @return Returns the object itself
#'
#' @keywords internal
print.mvoslr_power_object <- function(power_obj){

  if(length(power_obj$variable_parameter) == 0 |
     is.na(power_obj$variable_parameter) |
     is.null(power_obj$variable_parameter)){
    cat("Empirical power achieved:", power_obj$power, "\n")
  } else {
    cat("Empirical power achieved for all parameter choices:\n")
    print(power_obj$power)
  }

  invisible(power_obj)

}

#' Summary method for an object of class "mvoslr_power_object"
#'
#' @param power_obj Object of class "mvoslr_power_object"
#'
#' @return Returns the object itself
#'
#' @keywords internal
summary.mvoslr_power_object <- function(power_obj){

  # Show (interim) test decisions here
  cat("Results of multivariate one-sample log-rank test with varying ", attributes(power_obj)$variable_parameter, ".\n", sep = "")
  cat("Parameter assumes the following values:\n")
  print(attributes(power_obj)$parameter_values)

  if(attributes(power_obj)$variable_parameter == "a" & attributes(power_obj)$fixed_follow_up){
    cat("\nDuration of follow-up period was held fixed.\n")
  } else if(attributes(power_obj)$variable_parameter == "a" & attributes(power_obj)$fixed_follow_up){
    cat("\nOverall duration of trial was held fixed.\n")
  }
  cat("\n")

  if(length(power_obj$variable_parameter) == 0 |
     is.na(power_obj$variable_parameter) |
     is.null(power_obj$variable_parameter)){
    cat("Empirical power achieved after ",
        attributes(power_obj)$simulation_runs, " simulation runs: ",
        power_obj$power, "\n", sep = "")
  } else {
    cat("Empirical power achieved after ", attributes(power_obj)$simulation_runs,
        " simulation runs for all parameter choices:\n")
    print(power_obj$power)
  }

  cat("Rejection rates by stages:\n")
  print(power_obj$rejection_stages)

  # Show means if all parameters are held fixed
  if(is.matrix(power_obj$means)){
    cat("Means of multivariate process at all analysis dates:\n")
    print(power_obj$means)
  }

  invisible(power_obj)

}

#' Aggregates the results of several simulations. Before the aggregation itself, it is checked whether the objects are compatible in terms of the varied parameter and its values.
#'
#' @param power_obj_list List of objects of class "mvoslr_power_object"
#'
#' @return Object of type "mvoslr_power_object" in which the information of the objects in the list was aggregated
#'
#' @keywords internal
aggregate_mvoslr_power_object <- function(power_obj_list){

  # Check if information about power objects in the list is the same
  parameter_name_list <- lapply(power_obj_list, FUN = function(x) attributes(x)$variable_parameter)
  parameter_values_list <- lapply(power_obj_list, FUN = function(x) attributes(x)$parameter_values)
  follow_up_fixed_values <- lapply(power_obj_list, FUN = function(x) attributes(x)$follow_up_fixed)

  if( !(do.call(all.equal, parameter_name_list)) ){
    stop("Results of power calculations can't be aggregated due to incompatible parameter names.")
  }
  if( !(do.call(all.equal, parameter_values_list)) ){
    stop("Results of power calculations can't be aggregated due to incompatible parameter values.")
  }
  if( !(do.call(all.equal, follow_up_fixed_values)) ){
    stop("Results of power calculations can't be aggregated due to incompatible information about follow-up period.")
  }

  # Save attribute values of first list entry
  variable_parameter <- attributes(power_obj_list[[1]])$variable_parameter
  parameter_values <-attributes(power_obj_list[[1]])$parameter_values
  follow_up_fixed <- attributes(power_obj_list[[1]])$fixed_follow_up

  # Compute overall number of simulations and resulting weights
  runs_list <- lapply(power_obj_list, FUN = function(x) x$runs)
  overall_runs <- do.call(sum, runs_list)
  weights_list <- lapply(runs_list, FUN = function(n) n/overall_runs)

  # Collect information of each power object in lists
  power_list <- lapply(power_obj_list, FUN = function(x) x$power)
  rejection_stages_list <- lapply(power_obj_list, FUN = function(x) x$rejection_stages)
  means_list <- lapply(power_obj_list, FUN = function(x) x$means)
  variances_list <- lapply(power_obj_list, FUN = function(x) x$variances)

  # Multiply objects by the weight of the object
  weighted_power_list <- Map("*", power_list, weights_list)
  weighted_rejection_stages_list <- Map("*", rejection_stages_list, weights_list)
  weighted_means_list <- Map("*", means_list, weights_list)
  weighted_variances_list <- Map("*", variances_list, weights_list)

  # Sum up weighted objects to obtain means of the corresponding entities
  overall_power <- Reduce("+", weighted_power_list)
  overall_rejection_stages <- Reduce("+", weighted_rejection_stages_list)
  overall_means <- Reduce("+", weighted_means_list)
  overall_variances <- Reduce("+", weighted_variances_list)

  # Save averaged entities in new object
  overall_obj <- new_mvoslr_power_object(power = overall_power,
                                         rejection_stages = overall_rejection_stages,
                                         means = overall_means,
                                         variances = overall_variances,
                                         runs = overall_runs,
                                         variable_parameter = variable_parameter,
                                         parameter_values = parameter_values,
                                         follow_up_fixed = follow_up_fixed)

  return(overall_obj)

}
