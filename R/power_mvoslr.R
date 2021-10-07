#' Estimation of power of group-sequential multivariate one-sample log-rank test via simulation
#'
#' @param transition_matrix Matrix of transitions between states as in mstate package
#' @param model_type Reference multi-state model is either Markov (\code{model_type = "M"}) or Semi-Markov (\code{model_type = "SM"})
#' @param events List of (composite) events that shall be investigated
#' @param cum_hazard_functions_h0 Cumulative hazard functions for transitions in this model under null hypothesis
#' @param analysis_dates Vector of calendar dates of analyses
#' @param accrual_duration Duration of accrual period
#' @param sample_size Sample size for which power shall be estimated
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
#' @return List of 4:
#' \itemize{
#'   \item power - Overall power of the procedure under the specified alternative
#'   \item rejection_stages - Summary of stages in which the null hypothesis is rejected
#'   \item means - Mean value of the multivariate process from which test statistics are computed for all stages
#'   \item variances - Mean value of (co)variance estimates for the multivariate process for all stages
#' }
#'
#' @export
#'
#' @examples
#' #Setup of reference multi-state model (here: simple illness-death model)
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' number_of_trans_example <- dim(to.trans2(tmat_example))[1]
#' model_type_example <- "SM"
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' analysis_dates_example <- c(1, 2)
#' events_example <- list(c(2,3), c(3))
#' names(events_example) <- c("PFS", "OS")
#' accrual_duration_example <- 1
#' sample_size_example <- 100
#' #In this example, the alternative is specified via separate hazard ratios for each transition
#' hazard_ratios_example <- c(1.4, 1.2, 1.35)
#' power_mvoslr(transition_matrix = tmat_example, model_type = model_type_example,
#'              events = events_example, cum_hazard_functions_h0 = cum_hazards_example,
#'              analysis_dates = analysis_dates_example,
#'              accrual_duration = accrual_duration_example, sample_size = sample_size_example,
#'              hazard_ratios = hazard_ratios_example, simulation_runs = 10)
power_mvoslr <- function(transition_matrix, model_type, events, cum_hazard_functions_h0, analysis_dates, accrual_duration, sample_size,
                         hazard_ratios = NULL, cum_hazard_functions_alternative = NULL,
                         norm = "l2", boundaries = "obf", alpha = 0.05, weights = NULL, time_steps = 100,
                         simulation_runs = 1000){

  transitions <- mstate::to.trans2(transition_matrix)

  num_events <- length(events)
  num_transitions <- dim(transitions)[1]
  num_analyses <- length(analysis_dates)

  final_analysis <- max(analysis_dates)

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
  time <- seq(0, final_analysis, final_analysis/time_steps)

  cumhaz_alternative <- data.frame(time = rep(time, 3),
                                   Haz = unlist(lapply(cum_hazard_functions_alternative, function(x) x(time))),
                                   trans = rep(1:num_transitions, each = time_steps + 1))

  # Collect stagewise p-values, overall decisions, rejection stages and (raw) stagewise test-statistics
  p_collection <- matrix(NA, ncol = num_analyses, nrow = simulation_runs)
  decision_collection <- rep(NA, simulation_runs)
  rejection_stage_collection <- rep(NA, simulation_runs)
  stagewise_test_stat_collection <- array(NA, dim = c(num_events, num_analyses, simulation_runs))

  for(i in 1:simulation_runs){

    # simulation might fail for unknown reasons, repeat simulation function until it does not fail
    successful_simulation <- FALSE
    while(!successful_simulation){
      if(model_type == "M"){
        sim_data <- try(mstate::mssample(Haz = cumhaz_alternative,
                                         trans = transition_matrix,
                                         M = sample_size,
                                         output = "data"), silent = TRUE)
      } else if(model_type == "SM"){
        sim_data <- try(mstate::mssample(Haz = cumhaz_alternative,
                                         trans = transition_matrix,
                                         M = sample_size,
                                         clock = "reset",
                                         output = "data"), silent = TRUE)
      }
      if (class(sim_data) != "try-error") successful_simulation <- TRUE
    }

    ### Step 6: Transform simulation results

    sim_frame <- msdata_to_df(sim_data)

    ### Step 7: Simulate recruitment dates, append to simulation results, adapt observations according to censoring

    recruitment_dates <- runif(n = sample_size,
                               min = 0,
                               max = accrual_duration)

    ids_occurences <- table(sim_frame$id)
    sim_frame$recruitment_date <- rep(recruitment_dates,
                                      ids_occurences)
    sim_frame$censoring_date <- final_analysis - sim_frame$recruitment_date
    sim_frame$status <- sim_frame$status * (sim_frame$Tstop <= sim_frame$censoring_date)

    #adapt end and duration of observation period according to censoring
    sim_frame$Tstop <- pmin(sim_frame$Tstop, sim_frame$censoring_date)
    sim_frame$duration <- sim_frame$Tstop - sim_frame$Tstart

    sim_frame$censoring_date <- NULL

    result <- execution_mvoslr(msm_data = sim_frame, analysis_dates = analysis_dates, current_analysis = num_analyses,
                               transition_matrix = transition_matrix, cum_hazard_functions = cum_hazard_functions_h0,
                               model_type = model_type, events = events, norm = norm,
                               boundaries = boundaries, alpha = alpha, weights = weights)

    p_collection[i, ] <- result$stagewise_p_values
    decision_collection[i] <- ifelse(result$decision == "Accept H0", 0, 1)
    rejection_stage_collection[i] <- result$rejection_stage
    stagewise_test_stat_collection[ , , i] <- result$raw_martingale

  }

  power <- mean(decision_collection)
  rejection_stages <- table(rejection_stage_collection, exclude = NULL)/simulation_runs

  mean_summary <- apply(stagewise_test_stat_collection, MARGIN = c(1,2), FUN = mean)
  variance_summary <- array(NA, dim = c(num_events, num_events, num_analyses))
  for(a in 1:num_analyses){
    for(e1 in 1:num_events){
      for(e2 in e1:num_events){
        variance_summary[e1, e2, a] <-
          variance_summary[e2, e1, a] <- cov(stagewise_test_stat_collection[e1, a, ],
                                             stagewise_test_stat_collection[e2, a, ])
      }
    }
  }

  power_analysis <- list(power = power,
                         rejection_stages = rejection_stages,
                         means = mean_summary,
                         variances = variance_summary)

  return(power_analysis)
}

