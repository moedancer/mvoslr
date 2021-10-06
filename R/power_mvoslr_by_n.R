#' Estimation of power of group-sequential multivariate one-sample log-rank test via simulation with varying sample size
#'
#' @param transition_matrix Matrix of transitions between states as in mstate package
#' @param model_type Reference multi-state model is either Markov (\code{model_type = "M"}) or Semi-Markov (\code{model_type = "SM"})
#' @param events List of (composite) events that shall be investigated
#' @param cum_hazard_functions_h0 Cumulative hazard functions for transitions in this model under null hypothesis
#' @param analysis_dates Vector of calendar dates of analyses
#' @param accrual_duration Duration of accrual period
#' @param sample_sizes Vector of sample size for which power shall be estimated. As the accrual duration is kept fixed, the
#'                     recruitment speed differs for each sample size considered here!
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
#' @return List of 3:
#' \itemize{
#'   \item power - Overall power of the procedure under the specified alternative for each sample size
#'   \item rejection_stages - Summary of stages in which the null hypothesis is rejected for each sample size
#'   \item means - Mean value of the multivariate process from which test statistics are computed at each stage and for each sample size
#' }
#' For the sake of clarity of the output, mean estimated (co)variance matrices are not reported here.
#'
#' @export
#'
#' @examples
#' #Setup of reference multi-state model (here: simple illness-death model)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' number_of_trans_example <- dim(to.trans2(tmat_example))[1]
#' model_type_example <- "SM"
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' analysis_dates_example <- c(1, 2)
#' events_example <- list(c(2,3), c(3))
#' accrual_duration_example <- 1
#' sample_sizes_example <- seq(50, 150, 5)
#' #In this example, the alternative is specified via separate hazard ratios for each transition
#' hazard_ratios_example <- c(1.4, 1.2, 1.35)
#' power_mvoslr(transition_matrix = tmat_example, model_type = model_type_example,
#'              cum_hazard_functions_h0 = cum_hazards_example, analysis_dates = analysis_dates_example,
#'              accrual_duration = accrual_duration_example, sample_sizes = sample_sizes_example,
#'              hazard_ratios = hazard_ratios_example)
power_mvoslr_by_n <- function(transition_matrix, model_type, events, cum_hazard_functions_h0, analysis_dates,
                              accrual_duration, sample_sizes, hazard_ratios = NULL, cum_hazard_functions_alternative = NULL,
                              norm = "l2", boundaries = "obf", alpha = 0.05, weights = NULL, time_steps = 100,
                              simulation_runs = 1000){

  transitions <- to.trans2(transition_matrix)

  num_events <- length(events)
  num_transitions <- dim(transitions)[1]
  num_analyses <- length(analysis_dates)
  num_sample_sizes <- length(sample_sizes)

  final_analysis <- max(analysis_dates)

  max_sample_size <- max(sample_sizes)

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

  # Prepare arrays and matrices to collect stagewise p-values, overall decisions, rejection stages and (raw) stagewise test-statistics
  p_collection <- array(NA, dim = c(num_analyses, num_sample_sizes, simulation_runs))
  decision_collection <- matrix(NA, nrow = simulation_runs, ncol = num_sample_sizes)
  rejection_stage_collection <- matrix(NA, nrow = simulation_runs, ncol = num_sample_sizes)
  stagewise_test_stat_collection <- array(NA, dim = c(num_events, num_analyses, num_sample_sizes, simulation_runs))

  # Prepare names for aggregation of simulation results
  colnames(decision_collection) <- colnames(rejection_stage_collection) <- paste("n=", sample_sizes, sep = "")

  for(i in 1:simulation_runs){

    # Simulate trial data according to model under H0
    # Note: Simulation might fail for unknown reasons, repeat simulation function until it does not fail
    successful_simulation <- FALSE
    while(!successful_simulation){
      if(model_type == "M"){
        sim_data <- try(mssample(Haz = cumhaz_alternative,
                                 trans = transition_matrix,
                                 M = max_sample_size,
                                 output = "data"), silent = TRUE)
      } else if(model_type == "SM"){
        sim_data <- try(mssample(Haz = cumhaz_alternative,
                                 trans = transition_matrix,
                                 M = max_sample_size,
                                 clock = "reset",
                                 output = "data"), silent = TRUE)
      }
      if (class(sim_data) != "try-error") successful_simulation <- TRUE
    }

    # Transform simulation results

    sim_frame <- msdata_to_df(sim_data)

    # Simulate recruitment dates, append to simulation results, adapt observations according to censoring

    recruitment_dates <- runif(n = max_sample_size,
                               min = 0,
                               max = accrual_duration)

    ids_occurences <- table(sim_frame$id)
    sim_frame$recruitment_date <- rep(recruitment_dates,
                                      ids_occurences)
    sim_frame$censoring_date <- final_analysis - sim_frame$recruitment_date
    sim_frame$status <- sim_frame$status * (sim_frame$Tstop <= sim_frame$censoring_date)

    # Adapt end and duration of observation period according to censoring
    sim_frame$Tstop <- pmin(sim_frame$Tstop, sim_frame$censoring_date)
    sim_frame$duration <- sim_frame$Tstop - sim_frame$Tstart

    sim_frame$censoring_date <- NULL

    # Collect results of simulated trial
    result <- execution_mvoslr_by_n(msm_data = sim_frame, analysis_dates = analysis_dates, current_analysis = num_analyses,
                                    transition_matrix = transition_matrix, cum_hazard_functions = cum_hazard_functions_h0,
                                    model_type = model_type, events = events, sample_sizes = sample_sizes,
                                    norm = norm, boundaries = boundaries, alpha = alpha, weights = weights)

    p_collection[ , , i] <- result$stagewise_p_values
    decision_collection[i, ] <- ifelse(result$decision == "Accept H0", 0, 1)
    rejection_stage_collection[i, ] <- result$rejection_stage
    stagewise_test_stat_collection[ , , , i] <- result$raw_martingale

  }

  # calculate means over simulation runs
  power <- apply(decision_collection, 2, mean)
  rejection_stages <- apply(rejection_stage_collection, 2, function(stages){table(stages, exclude = NULL)/simulation_runs})

  # Helper function to combine overview over rejection stages for different sample sizes into one matrix
  custom_rbind <- function(data_a, data_b){
    if(!is.matrix(data_a)) data_a <- t(as.matrix(data_a))
    if(!is.matrix(data_b)) data_b <- t(as.matrix(data_b))
    colnames(data_a) <- replace(colnames(data_a), which(is.na(colnames(data_a))), "Acceptance")
    colnames(data_b) <- replace(colnames(data_b), which(is.na(colnames(data_b))), "Acceptance")
    diff_1 <- setdiff(colnames(data_a), colnames(data_b))
    diff_2 <- setdiff(colnames(data_b), colnames(data_a))
    for(x in diff_1){
      data_b <- cbind(data_b, 0)
      colnames(data_b)[dim(data_b)[2]] <- x
    }
    for(x in diff_2){
      data_a <- cbind(data_a, 0)
      colnames(data_a)[dim(data_a)[2]] <- x
    }
    data_b <- data_b[,colnames(data_a)]
    return(rbind(data_a, data_b))
  }

  # Create matrix for overview over rejection stages (if neccessary)
  if(is.list(rejection_stages)){
    rejection_stages <- Reduce(custom_rbind, rejection_stages)
    rownames(rejection_stages) <- paste("n=", sample_sizes, sep = "")
    rejection_stages <- t(rejection_stages)
    rejection_stages <- rejection_stages[sort(rownames(rejection_stages)), ]
  } else {
    rownames(rejection_stages) <- replace(rownames(rejection_stages), which(is.na(rownames(rejection_stages))), "Acceptance")
  }

  # Compute means of raw martingales for each analysis and choice of sample size
  mean_summary <- apply(stagewise_test_stat_collection, MARGIN = c(1,2,3), FUN = mean)
  if(!is.null(names(events))){
    dimnames(mean_summary)[[1]] <- names(events)
  } else {
    dimnames(mean_summary)[[1]] <- paste("Event", 1:num_events, sep = " ")
  }
  dimnames(mean_summary)[[2]] <- paste("Analysis", 1:num_analyses, sep = " ")
  dimnames(mean_summary)[[3]] <- paste("n=", sample_sizes, sep = "")

  # Overview of estimated covariance matrices is omitted (would result in 4-dimensional array)

  power_analysis <- list(power = power,
                         rejection_stages = rejection_stages,
                         means = mean_summary)

  return(power_analysis)
}
