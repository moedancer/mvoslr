#' Estimation of power of group-sequential multivariate one-sample log-rank test via simulation with varying accrual duration
#'
#' In this case the calendar dates of interim analyses are fixed and the minimum follow-up length are fixed.
#' As the accrual duration can be varied, the calendar date of the final analysis varies.
#'
#' @param reference_model Specification of the reference model against which the new data is tested. Should be an object of class "reference_model"
#' @param events List of (composite) events that shall be investigated
#' @param interim_analysis_dates Vector of calendar dates of interim analyses. Calendar date of final analysis depends on the respective accrual duration
#' @param accrual_durations Vector of durations of accrual period
#' @param follow_up Length of minimum follow-up
#' @param recruitment_speed Number of patients that can be recruited during one unit of time.
#' @param hazard_ratios Specification of alternative hypothesis. One can either choose a single hazard ratio which then applies
#'                      to each transition or a single hazard ratio for each transition in the model.
#' @param cum_hazard_functions_alternative If non-proportional transition intensities are anticipated, the transition intensities
#'                                         for all transitions, need to be specified here. If this argument is specified, it
#'                                         overrules the argument \code{hazard_ratios}.
#' @param norm Use either \eqn{L^2}-norm (\code{norm = "l2"}, default value) or \eqn{L^\infty}-norm (\code{norm = "linf"}) of vector of test statistics to compute stagewise p-values
#' @param boundaries Use either O'Brien-Fleming'S (\code{boundaries = "obf"}, default value) or Pocock's (\code{boundaries = "pocock"}) sequential decision boundaries
#' @param alpha Choose type I error rate (default value = 0.05)
#' @param weights Choose weights for inverse normal combination of stagewise p-values. Sum of squared values needs to sum up to 1.
#' @param time_steps As the multi-state model data will be simulated with the \code{mstate}-package, the transition intensities need
#'                   to be discretized. The number of time steps in which the intensities will be discretized can be specified here.
#'                   The time horizon over which the discretization happens is the calndar date of the last analysis.
#' @param simulation_runs The number of simulation runs to estimate the power can be specified here.
#'
#' @return Object of class "mvoslr_power_object"
#' For the sake of clarity of the output, mean estimated (co)variance matrices are not reported here.
#'
#' @export
#'
#' @examples
#' #Setup of reference multi-state model (here: simple illness-death model)
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' number_of_trans_example <- dim(to.trans2(tmat_example))[1]
#' model_type_example <- "SM"
#' events_example <- list(c(2,3), c(3))
#' names(events_example) <- c("PFS", "OS")
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' reference_model_example <- new_reference_model(transition_matrix = tmat_example,
#'                                                intensities = cum_hazards_example,
#'                                                type = model_type_example)
#' interim_analysis_dates_example <- 1
#' accrual_durations_example <- seq(0.5, 1.5, 0.2)
#' follow_up_example <- 1
#' recruitment_speed_example <- 100
#' #In this example, the alternative is specified via separate hazard ratios for each transition
#' hazard_ratios_example <- c(1.4, 1.2, 1.35)
#' power_mvoslr_fixed_fu(reference_model = reference_model_example,
#'                       events = events_example,
#'                       interim_analysis_dates = interim_analysis_dates_example,
#'                       accrual_durations = accrual_durations_example,
#'                       follow_up = follow_up_example,
#'                       recruitment_speed = recruitment_speed_example,
#'                       hazard_ratios = hazard_ratios_example, simulation_runs = 10)
power_mvoslr_fixed_fu <- function(reference_model, events, interim_analysis_dates,
                                  accrual_durations, follow_up, recruitment_speed, hazard_ratios = NULL, cum_hazard_functions_alternative = NULL,
                                  norm = "l2", boundaries = "obf", alpha = 0.05, weights = NULL, time_steps = 100,
                                  simulation_runs = 1000){

  # Unpack information from reference model
  transition_matrix <- reference_model$transition_matrix
  cum_hazard_functions_h0 <- reference_model$intensities
  model_type <- attributes(reference_model)$type

  transitions <- mstate::to.trans2(transition_matrix)

  num_events <- length(events)
  num_transitions <- dim(transitions)[1]
  num_analyses <- length(interim_analysis_dates) + 1
  num_accrual_durations <- length(accrual_durations)

  longest_accrual <- max(accrual_durations)
  final_analysis <- longest_accrual + follow_up

  sample_size <- longest_accrual * recruitment_speed

  # Check if alternative is specified
  # If it is specified in terms of hazard rations, compute corresponding cumulative hazard functions
  if(is.null(hazard_ratios) & is.null(cum_hazard_functions_alternative)){
    return("Specify an alternative! Either in terms of hazard ratios or by supplying separate cumulative hazard functions.")
  } else if(!is.null(hazard_ratios) & length(hazard_ratios) == 1){
    cum_hazard_functions_alternative <- lapply(cum_hazard_functions_h0,
                                               function(f){
                                                 function(x) (1/hazard_ratios) * do.call(f, list(x))
                                               })
  } else if(length(hazard_ratios) > 1){
    cum_hazard_functions_alternative <- vector(mode = "list", length = num_transitions)
    for(trans in 1:num_transitions){
      cum_hazard_functions_alternative[[trans]] <- function(x){
        (1/hazard_ratios[trans]) * do.call(cum_hazard_functions_h0[[trans]], list(x))
      }
    }
  }

  # Discretize hazards under planning alternative for simulation with mstate
  cumhaz_alternative <- discretize_functions(cum_hazard_functions_alternative,
                                             final_analysis, time_steps = time_steps)

  # Prepare arrays and matrices to collect stagewise p-values, overall decisions, rejection stages and (raw) stagewise test-statistics
  p_collection <- array(NA, dim = c(num_analyses, num_accrual_durations, simulation_runs))
  decision_collection <- matrix(NA, nrow = simulation_runs, ncol = num_accrual_durations)
  rejection_stage_collection <- matrix(NA, nrow = simulation_runs, ncol = num_accrual_durations)
  stagewise_test_stat_collection <- array(NA, dim = c(num_events, num_analyses, num_accrual_durations, simulation_runs))

  # Prepare names for aggregation of simulation results
  colnames(decision_collection) <- colnames(rejection_stage_collection) <- paste("a=", accrual_durations, sep = "")

  for(i in 1:simulation_runs){

    sim_frame <- simulate_msm(transition_matrix, model_type, cumhaz_alternative, sample_size)

    sim_frame <- msm_to_trial_data(sim_frame, longest_accrual,
                                   follow_up)

    # Collect results of simulated trial
    result <- execution_mvoslr_fixed_fu(msm_data = sim_frame, interim_analysis_dates = interim_analysis_dates,
                                        current_analysis = num_analyses, follow_up = follow_up,
                                        reference_model = reference_model, events = events, accrual_durations = accrual_durations,
                                        norm = norm, boundaries = boundaries, alpha = alpha, weights = weights)

    p_collection[ , , i] <- result$stagewise_p_values
    decision_collection[i, ] <- ifelse(is.na(result$rejection_stage), 0, 1)
    rejection_stage_collection[i, ] <- result$rejection_stage
    stagewise_test_stat_collection[ , , , i] <- result$raw_process

  }

  # calculate means over simulation runs
  power <- apply(decision_collection, 2, mean)
  rejection_stage_collection[is.na(rejection_stage_collection)] <- 0
  rejection_stages <- apply(rejection_stage_collection, 2,
                            function(stages){table(factor(stages, levels = 0:num_analyses), exclude = NULL)/simulation_runs})

  rownames(rejection_stages) <- replace(rownames(rejection_stages), which(rownames(rejection_stages) == 0), "Acceptance")

  # Compute means of raw martingales for each analysis and choice of accrual period
  mean_summary <- apply(stagewise_test_stat_collection, MARGIN = c(1,2,3), FUN = mean)

  # Overview of estimated covariance matrices is omitted (would result in 4-dimensional array)

  output_obj <- new_mvoslr_power_object(power = power,
                                        rejection_stages = rejection_stages,
                                        means = mean_summary,
                                        variable_parameter = "a",
                                        parameter_values = accrual_durations,
                                        runs = simulation_runs,
                                        follow_up_fixed = TRUE)

  return(output_obj)
}


#' Estimation of power of group-sequential multivariate one-sample log-rank test via simulation with varying accrual duration
#'
#' In this case the calendar dates of interim analyses are fixed and the minimum follow-up length are fixed.
#' As the accrual duration can be varied, the calendar date of the final analysis varies.
#'
#' @param reference_model Specification of the reference model against which the new data is tested. Should be an object of class "reference_model"
#' @param events List of (composite) events that shall be investigated
#' @param interim_analysis_dates Vector of calendar dates of interim analyses. Calendar date of final analysis depends on the respective accrual duration
#' @param accrual_durations Vector of durations of accrual period
#' @param follow_up Length of minimum follow-up
#' @param recruitment_speed Number of patients that can be recruited during one unit of time.
#' @param hazard_ratios Specification of alternative hypothesis. One can either choose a single hazard ratio which then applies
#'                      to each transition or a single hazard ratio for each transition in the model.
#' @param cum_hazard_functions_alternative If non-proportional transition intensities are anticipated, the transition intensities
#'                                         for all transitions, need to be specified here. If this argument is specified, it
#'                                         overrules the argument \code{hazard_ratios}.
#' @param norm Use either \eqn{L^2}-norm (\code{norm = "l2"}, default value) or \eqn{L^\infty}-norm (\code{norm = "linf"}) of vector of test statistics to compute stagewise p-values
#' @param boundaries Use either O'Brien-Fleming'S (\code{boundaries = "obf"}, default value) or Pocock's (\code{boundaries = "pocock"}) sequential decision boundaries
#' @param alpha Choose type I error rate (default value = 0.05)
#' @param weights Choose weights for inverse normal combination of stagewise p-values. Sum of squared values needs to sum up to 1.
#' @param time_steps As the multi-state model data will be simulated with the \code{mstate}-package, the transition intensities need
#'                   to be discretized. The number of time steps in which the intensities will be discretized can be specified here.
#'                   The time horizon over which the discretization happens is the calndar date of the last analysis.
#' @param simulation_runs The number of simulation runs to estimate the power can be specified here.
#' @param cores Number of cores to be used for the simulation
#'
#' @return Object of class "mvoslr_power_object"
#' For the sake of clarity of the output, mean estimated (co)variance matrices are not reported here.
#'
#' @export
#'
#' @examples
#' #Setup of reference multi-state model (here: simple illness-death model)
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' number_of_trans_example <- dim(to.trans2(tmat_example))[1]
#' model_type_example <- "SM"
#' events_example <- list(c(2,3), c(3))
#' names(events_example) <- c("PFS", "OS")
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' reference_model_example <- new_reference_model(transition_matrix = tmat_example,
#'                                                intensities = cum_hazards_example,
#'                                                type = model_type_example)
#' interim_analysis_dates_example <- 1
#' accrual_durations_example <- seq(0.5, 1.5, 0.2)
#' follow_up_example <- 1
#' recruitment_speed_example <- 100
#' #In this example, the alternative is specified via separate hazard ratios for each transition
#' hazard_ratios_example <- c(1.4, 1.2, 1.35)
#' power_mvoslr_fixed_fu_par(reference_model = reference_model_example,
#'                           events = events_example,
#'                           interim_analysis_dates = interim_analysis_dates_example,
#'                           accrual_durations = accrual_durations_example,
#'                           follow_up = follow_up_example,
#'                           recruitment_speed = recruitment_speed_example,
#'                           hazard_ratios = hazard_ratios_example, simulation_runs = 10, cores = 2)
power_mvoslr_fixed_fu_par <- function(reference_model, events, interim_analysis_dates,
                                      accrual_durations, follow_up, recruitment_speed, hazard_ratios = NULL, cum_hazard_functions_alternative = NULL,
                                      norm = "l2", boundaries = "obf", alpha = 0.05, weights = NULL, time_steps = 100,
                                      simulation_runs = 1000, cores = NULL){
  if(!requireNamespace("foreach", quietly = TRUE) |
     !requireNamespace("parallel", quietly = TRUE) |
     !requireNamespace("doParallel", quietly = TRUE)){

    warning("Parallelised version cannot be used as suggested packages aren't installed.
            Standard power calculation function is used instead.")

    result <- power_mvoslr_fixed_fu(reference_model = reference_model, events = events,
                                    interim_analysis_dates = interim_analysis_dates,
                                    accrual_durations = accrual_durations,
                                    follow_up = follow_up, recruitment_speed = recruitment_speed,
                                    hazard_ratios = hazard_ratios,
                                    cum_hazard_functions_alternative = cum_hazard_functions_alternative,
                                    norm = norm, boundaries = boundaries, alpha = alpha, weights = weights,
                                    time_steps = time_steps, simulation_runs = simulation_runs)

  } else {

    if(is.null(cores)) cores <- parallel::detectCores()

    # Distribute number of simulations to cores
    distributed_runs <- rep(floor(simulation_runs/cores), cores)
    rest <- simulation_runs %%cores
    distributed_runs[1:rest] <- distributed_runs[1:rest] + 1

    doParallel::registerDoParallel(cores)

    # Define counter to pass checks
    i <- NULL

    results_list <- foreach::"%dopar%"(foreach::foreach(i = 1:cores, .combine = list),

                                      power_mvoslr_fixed_fu(reference_model = reference_model, events = events,
                                                            interim_analysis_dates = interim_analysis_dates,
                                                            accrual_durations = accrual_durations,
                                                            follow_up = follow_up, recruitment_speed = recruitment_speed,
                                                            hazard_ratios = hazard_ratios,
                                                            cum_hazard_functions_alternative = cum_hazard_functions_alternative,
                                                            norm = norm, boundaries = boundaries, alpha = alpha, weights = weights,
                                                            time_steps = time_steps, simulation_runs = distributed_runs[i])

    )

    result <- aggregate_mvoslr_power_object(results_list)

  }

  return(result)

}



