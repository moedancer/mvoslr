#' Simulation of trial data given model specifications and trial design
#'
#' @param transition_matrix matrix of transitions between states as in mstate package
#' @param model_type Multi-state model is either Markov (\code{model_type = "M"}) or Semi-Markov (\code{model_type = "SM"})
#' @param cum_hazard_functions Cumulative hazard functions for transitions in this model
#' @param accrual_duration Duration of accrual period
#' @param follow_up_duration Duration of follow-up period
#' @param sample_size Number of patients to simulate
#' @param time_steps As the multi-state model data will be simulated with the \code{mstate}-package, the transition intensities need
#'                   to be discretized. The number of time steps in which the intensities will be discretized can be specified here.
#'                   The time horizon over which the discretization happens is the calndar date of the last analysis.
#'
#' @return Data frame with simulated trial data
#'
#' @export
#'
#' @examples
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' model_type_example <- "SM"
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' accrual_duration_example <- 1
#' follow_up_duration_example <- 1
#' sample_size_example <- 100
#' simulate_trial(transition_matrix = tmat_example, model_type = model_type_example,
#'                cum_hazard_functions = cum_hazards_example,
#'                accrual_duration = accrual_duration_example,
#'                follow_up_duration = follow_up_duration_example,
#'                sample_size = sample_size_example,)
simulate_trial <- function(transition_matrix, model_type, cum_hazard_functions,
                           accrual_duration, follow_up_duration, sample_size,
                           time_steps = 100){

  transitions <- mstate::to.trans2(transition_matrix)

  cumhaz_frame <- discretize_functions(function_list = cum_hazard_functions,
                                       max_argument = accrual_duration + follow_up_duration,
                                       time_steps = time_steps)

  # Simulate multi-state data with mstate
  # Simulation might fail for unknown reasons, repeat simulation function until it does not fail
  sim_frame <- simulate_msm(transition_matrix, model_type, cumhaz_frame, sample_size)

  sim_frame <- msm_to_trial_data(sim_frame, accrual_duration, follow_up_duration)

  return(sim_frame)

}
