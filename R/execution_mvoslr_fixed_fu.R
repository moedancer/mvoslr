#' Execution of interim and final analyses of multivariate group-sequential one-sample log-rank test with varying accrual duration
#' and a fixed follow-up after end of the accrual period
#'
#' For a given data set, this function can efficiently compute the results of the group-sequential multivariate one-sample
#' log-rank test for different durations of the accrual duration. Calendar dates of interim_analyses are fixed and the calendar
#' date of the final analysis is set to the duration of the follow-up period plus the duration of the respective accrual period.
#' This function is needed to efficiently estimate the power of such trials for different sample sizes via simulation.
#'
#' @param msm_data Data frame containing information about any (possible) transition of patients in the study with the following columns:
#' \itemize{
#'   \item id - ID of the patient to whom this entry belongs
#'   \item Tstart - Start of the observation period for this transition
#'   \item Tstop - Stop of the observation period for this transition
#'   \item duration - Duration of the observation period for this transition
#'   \item from - Current state of the patient
#'   \item to - State to which a transition is possible
#'   \item status - Transition could (1) or could not be (0) observed
#'   \item trans - Number of this transition (according to \code{to.trans2(transition_matrix)})
#'   \item recruitment_date - Calendar date of recruitment of this patient
#' }
#' Maximum number of sample sizes under consideration should be at most the number of subjects in this data set
#' @param interim_analysis_dates Vector of calendar dates of interim analyses
#' @param follow_up Duration of the follow-up period (after stop of recruitment)
#' @param current_analysis Number of the current analysis
#' @param reference_model Specification of the reference model against which the new data is tested. Should be an object of class "reference_model"
#' @param events List of (composite) events that shall be investigated
#' @param accrual_durations Vector of durations of the accrual period that shall be investigated
#' @param norm Use either \eqn{L^2}-norm (\code{norm = "l2"}, default value) or \eqn{L^\infty}-norm (\code{norm = "linf"}) of vector of test statistics to compute stagewise p-values
#' @param boundaries Use either O'Brien-Fleming'S (\code{boundaries = "obf"}, default value) or Pocock's (\code{boundaries = "pocock"}) sequential decision boundaries
#' @param alpha Choose type I error rate (default value = 0.05)
#' @param weights Choose weights for inverse normal combination of stagewise p-values. Sum of squared values needs to sum up to 1.
#'
#' @return List of 4:
#' \itemize{
#'   \item raw_martingale - Values of the multivariate process at all stages for all sample sizes
#'   \item stagewise_p_values - Stagewise p-values for all sample sizes
#'   \item decision - Decision to be made at the current analysis for each sample size
#'   \item decision_stage - Number of analysis at which the null hypothesis is rejected (\code{NA} if not applicable)
#'                          for each sample size
#' }
#'
#' @import stats
#'
#' @keywords internal
#'
execution_mvoslr_fixed_fu <- function(msm_data, interim_analysis_dates, follow_up, current_analysis = NULL, reference_model,
                                      events, accrual_durations, norm = "l2", boundaries = "obf", alpha = 0.05, weights = NULL){

  # Unpack information from reference model
  transition_matrix <- reference_model$transition_matrix
  cum_hazard_functions <- reference_model$intensities
  model_type <- attributes(reference_model)$type

  # Check some arguments of the function
  if(!model_type %in% c("M", "SM")){
    cat("Parameter model needs to be one of \"M\" (Markov model) or \"SM\" (Semi-Markov model)!" )
  } else if(!norm %in% c("l2", "linf")) {
    cat("Chosen norm need to be either \"l2\" or \"linf\"!")
  } else if(!boundaries %in% c("obf", "pocock")) {
    cat("Sequential boundaries needs to be one of \"obf\" or \"pocock\"!")
  } else if(!is.null(weights) & abs(1 - sum(weights^2)) > 0.00001){
    cat("Squared weights do not sum up to one!")
  } else if(!is.null(weights) & length(weights)!=length(interim_analysis_dates) + 1){
    cat("Length of weight vector does not correspond to number of analyses!")
  } else if(min(accrual_durations) + follow_up < max(interim_analysis_dates)){
    cat("Last interim analysis is not included in overall time horizon for at least one choice of accrual duration")
  } else {

    number_accrual_durations <- length(accrual_durations)
    number_event_types       <- length(events)
    number_of_analyses       <- length(interim_analysis_dates) + 1 # Note: time of final analysis is not fixed here!

    longest_accrual <- max(accrual_durations)

    transitions <- mstate::to.trans2(transition_matrix)
    number_of_trans <- dim(transitions)[1]

    # Use equal weights if not specified otherwise
    if(is.null(weights)) weights <- rep(1/sqrt(number_of_analyses), number_of_analyses)

    # Current analysis is set to overall number of analyses if not specified otherwise
    if(is.null(current_analysis)) current_analysis <- number_of_analyses

    # Create arrays to store cumulative and incremental events and accumulated hazards
    events_count  <- array(0, c(number_event_types, number_event_types, current_analysis, number_accrual_durations))
    events_hazard <- array(0, c(number_event_types, number_event_types, current_analysis, number_accrual_durations))

    events_count_inc  <- array(0, c(number_event_types, number_event_types, current_analysis, number_accrual_durations))
    events_hazard_inc <- array(0, c(number_event_types, number_event_types, current_analysis, number_accrual_durations))

    # Create matrices/vectors for martingale, test statistics and p-values
    mv_martingale <- array(0, dim = c(number_event_types, current_analysis, number_accrual_durations))
    mv_martingale_inc <- array(0, dim = c(number_event_types, current_analysis, number_accrual_durations))
    mv_test_statistic <- array(0, dim = c(number_event_types, current_analysis, number_accrual_durations))
    uv_test_statistic <- matrix(0, nrow = current_analysis, ncol = number_accrual_durations)

    # Assign names to arrays for better readability of outputs
    if(is.null(names(events))) names(events) <- paste("Event", 1:number_event_types, sep = " ")

    dimnames(mv_martingale)[[1]] <- names(events)
    dimnames(mv_test_statistic)[[1]] <- names(events)
    dimnames(events_count)[[1]] <- dimnames(events_count)[[2]] <- names(events)
    dimnames(events_hazard)[[1]] <- dimnames(events_hazard)[[2]] <- names(events)
    dimnames(events_count_inc)[[1]] <- dimnames(events_count_inc)[[2]] <- names(events)
    dimnames(events_hazard_inc)[[1]] <- dimnames(events_hazard_inc)[[2]] <- names(events)

    # Compute test statistics for each analysis date
    # Note: Date of last analysis isn't fixed and needs special treatment

    for(number_of_analysis in 1:(min(current_analysis, number_of_analyses - 1))){

      analysis_date <- interim_analysis_dates[number_of_analysis]

      # Introduce "..._temp" data.frame as each analysis implies different censoring pattern
      msm_data_temp <- msm_to_trial_data(msm_data, min(analysis_date, longest_accrual),
                                         ifelse(longest_accrual > analysis_date,
                                                0, analysis_date - longest_accrual))

      # Calculate accumulated hazards for each transition
      # Create column to enter accumulated hazards
      msm_data_temp$acc_haz <- rep(0, dim(msm_data_temp)[1])

      for(i in 1:number_of_trans){

        # Find correct hazard function
        current_cumhaz <- cum_hazard_functions[[i]]

        current_start <- transitions$from[i]
        current_end <- transitions$to[i]

        # Find entries belonging to the current transition
        current_transitions <- msm_data_temp[which(msm_data_temp$from == current_start &
                                                     msm_data_temp$to == current_end), ]

        # Compute accumulated hazard for this transition (according to the model type)
        if(model_type == "M"){
          current_acc_hazards <- current_cumhaz(current_transitions$Tstop) -
            current_cumhaz(current_transitions$Tstart)
        } else if(model_type == "SM"){
          current_acc_hazards <- current_cumhaz(current_transitions$duration)
        }
        msm_data_temp$acc_haz[which(msm_data_temp$from == current_start &
                                      msm_data_temp$to == current_end)] <- current_acc_hazards

      }

      # Compute value of counting process and overall accumulated hazard for each (composite) event (if i =j)
      #   and joint occurence of two events and overall accumulated hazard for this (if i != j)
      for(i in 1:number_event_types){

        for(j in i:number_event_types){

          event_1 <- events[[i]]
          event_2 <- events[[j]]
          union_event <- sort(unique(c(event_1, event_2)))
          intersect_event <- intersect(event_1, event_2)

          if(length(intersect_event) > 0){

            # Extract relevant information for current event
            relevant_obs <- msm_data_temp[which(msm_data_temp$to %in% union_event), ]
            true_events <- relevant_obs[which(relevant_obs$status == 1), ]

            # Compute event time
            event_times <- aggregate(Tstop ~ id, true_events, function(x) min(x))

            # Find IDs of patients without event
            all_patients <- unique(msm_data_temp$id)
            eventless_patients <- all_patients[which(!(all_patients %in% event_times$id))]
            number_eventless_patients <- length(eventless_patients)

            # Add dummy rows for those patients
            if(number_eventless_patients > 0){
              dummy_rows <- data.frame(id = eventless_patients, Tstop = rep(Inf, number_eventless_patients))
              event_times <- rbind(event_times, dummy_rows)
            }

            event_times <- event_times[order(event_times$id),]

            # Extract observations until occurence of event
            event_times_rep <-  event_times$Tstop[match(relevant_obs$id, event_times$id)]

            for(a_index in 1:number_accrual_durations){

              # Need "_temp" object for computation for different accrual durations
              # Note: we are interested in the history of the joint (intersection) event before union event
              relevant_obs_temp <- relevant_obs[which(relevant_obs$Tstop <= event_times_rep &
                                                        relevant_obs$to %in% intersect_event &
                                                        relevant_obs$recruitment_date <= accrual_durations[a_index]), ]

              # Count events and sum up hazards for (joint events)
              sum_intersect_events <- sum(relevant_obs_temp$status)
              sum_accumulated_hazard <- sum(relevant_obs_temp$acc_haz)

              events_count[i, j, number_of_analysis, a_index] <-
                events_count[j, i, number_of_analysis, a_index]  <- sum_intersect_events
              events_hazard[i, j, number_of_analysis, a_index] <-
                events_hazard[j, i, number_of_analysis, a_index] <- sum_accumulated_hazard

            }
          }
        }
      }

      # Compute value of M0 for each event
      mv_martingale[,number_of_analysis,] <-
        apply(events_count[,,number_of_analysis,], 3, diag) - apply(events_hazard[,,number_of_analysis,], 3, diag)

      # Compute increments of event count and cumulative hazard data and M0
      if(number_of_analysis == 1){

        events_count_inc[,,1,] <- events_count[,,1,]
        events_hazard_inc[,,1,] <- events_hazard[,,1,]
        mv_martingale_inc[,1,] <- mv_martingale[,1,]

      } else {

        events_count_inc[,,number_of_analysis,] <- events_count[,,number_of_analysis,] - events_count[,,number_of_analysis-1,]
        events_hazard_inc[,,number_of_analysis,] <- events_hazard[,,number_of_analysis,] - events_hazard[,,number_of_analysis-1,]
        mv_martingale_inc[,number_of_analysis,] <- mv_martingale[,number_of_analysis,] - mv_martingale[,number_of_analysis-1,]

      }

    }

    # Special treatment of last analysis as timing of last analysis depends on length of accrual duration
    if(current_analysis == number_of_analyses){

      for(a_index in 1:number_accrual_durations){

        analysis_date <- accrual_durations[a_index] + follow_up

        # Introduce "..._temp" data.frame as each analysis implies different censoring pattern
        msm_data_temp <- msm_to_trial_data(msm_data, accrual_durations[a_index],
                                           follow_up)

        # Calculate accumulated hazards for each transition
        # Create column to enter accumulated hazards
        msm_data_temp$acc_haz <- rep(0, dim(msm_data_temp)[1])

        for(i in 1:number_of_trans){

          # Find correct hazard function
          current_cumhaz <- cum_hazard_functions[[i]]

          current_start <- transitions$from[i]
          current_end <- transitions$to[i]

          # Find entries belonging to the current transition
          current_transitions <- msm_data_temp[which(msm_data_temp$from == current_start &
                                                       msm_data_temp$to == current_end), ]

          # Compute accumulated hazard for this transition (according to the model type)
          if(model_type == "M"){
            current_acc_hazards <- current_cumhaz(current_transitions$Tstop) -
              current_cumhaz(current_transitions$Tstart)
          } else if(model_type == "SM"){
            current_acc_hazards <- current_cumhaz(current_transitions$duration)
          }
          msm_data_temp$acc_haz[which(msm_data_temp$from == current_start &
                                        msm_data_temp$to == current_end)] <- current_acc_hazards

        }

        # Compute value of counting process and overall accumulated hazard for each (composite) event (if i =j)
        #   and joint occurence of two events and overall accumulated hazard for this (if i != j)
        for(i in 1:number_event_types){

          for(j in i:number_event_types){

            event_1 <- events[[i]]
            event_2 <- events[[j]]
            union_event <- sort(unique(c(event_1, event_2)))
            intersect_event <- intersect(event_1, event_2)

            if(length(intersect_event) > 0){

              # Extract relevant information for current event
              relevant_obs <- msm_data_temp[which(msm_data_temp$to %in% union_event), ]
              true_events <- relevant_obs[which(relevant_obs$status == 1), ]

              # Compute event time
              event_times <- aggregate(Tstop ~ id, true_events, function(x) min(x))

              # Find IDs of patients without event
              all_patients <- unique(msm_data_temp$id)
              eventless_patients <- all_patients[which(!(all_patients %in% event_times$id))]
              number_eventless_patients <- length(eventless_patients)

              # Add dummy rows for those patients
              if(number_eventless_patients > 0){
                dummy_rows <- data.frame(id = eventless_patients, Tstop = rep(Inf, number_eventless_patients))
                event_times <- rbind(event_times, dummy_rows)
              }

              event_times <- event_times[order(event_times$id),]

              # Extract observations until occurence of event
              event_times_rep <-  event_times$Tstop[match(relevant_obs$id, event_times$id)]

              # Need "_temp" object for computation for different accrual durations
              # Note: we are interested in the history of the joint (intersection) event before union event
              relevant_obs_temp <- relevant_obs[which(relevant_obs$Tstop <= event_times_rep &
                                                        relevant_obs$to %in% intersect_event &
                                                        relevant_obs$recruitment_date <= accrual_durations[a_index]), ]

              # Count events and sum up hazards for (joint events)
              sum_intersect_events <- sum(relevant_obs_temp$status)
              sum_accumulated_hazard <- sum(relevant_obs_temp$acc_haz)

              events_count[i, j, number_of_analyses, a_index] <-
                events_count[j, i, number_of_analyses, a_index]  <- sum_intersect_events
              events_hazard[i, j, number_of_analyses, a_index] <-
                events_hazard[j, i, number_of_analyses, a_index] <- sum_accumulated_hazard

            }
          }
        }
      }

      # Compute value of M0 for each event
      mv_martingale[,number_of_analyses,] <-
        apply(events_count[,,number_of_analyses,], 3, diag) - apply(events_hazard[,,number_of_analyses,], 3, diag)

      # Compute increments of event count and cumulative hazard data and M0
      if(number_of_analyses == 1){

        events_count_inc[,,1,] <- events_count[,,1,]
        events_hazard_inc[,,1,] <- events_hazard[,,1,]
        mv_martingale_inc[,1,] <- mv_martingale[,1,]

      } else {

        events_count_inc[,,number_of_analyses,] <- events_count[,,number_of_analyses,] - events_count[,,number_of_analyses-1,]
        events_hazard_inc[,,number_of_analyses,] <- events_hazard[,,number_of_analyses,] - events_hazard[,,number_of_analyses-1,]
        mv_martingale_inc[,number_of_analyses,] <- mv_martingale[,number_of_analyses,] - mv_martingale[,number_of_analyses-1,]

      }
    }

    # Estimate covariance matrix according to the (generalized) suggestion from Wu (2014)
    covariance_estimate <- 0.5 * events_count + 0.5 * events_hazard
    covariance_estimate_inc <- 0.5 * events_count_inc + 0.5 * events_hazard_inc

    # Compute Cholesky factorization to standardize multivariate process M0
    standardisation_matrices_raw <- apply(covariance_estimate_inc, c(3, 4), function(matrix) solve(t(chol(matrix))))
    standardisation_matrices <- array(standardisation_matrices_raw,
                                      dim = c(number_event_types, number_event_types, current_analysis, number_accrual_durations))

    mv_test_statistic <- array(0, dim=c(number_event_types, current_analysis, number_accrual_durations))
    for(i in 1:current_analysis){
      for(a_index in 1:number_accrual_durations){
        mv_test_statistic[, i, a_index] <- standardisation_matrices[, , i, a_index] %*% mv_martingale_inc[, i, a_index]
      }
    }

    if(!is.null(names(events))){
      dimnames(mv_test_statistic)[[1]] <- names(events)
    }

    # Compute univariate test statistics and stagewise p-values according to chosen norm (L2 or Linf)
    # Note: In "uv_test_statistic", stages are given as rows, sample sizes are given as columns
    if(norm == "l2"){
      uv_test_statistic <- apply(mv_test_statistic, c(2, 3), function(vec) sqrt(sum(vec^2)))
      stagewise_p_values <- 1 - pchisq(uv_test_statistic^2, df = number_event_types)
    } else if(norm == "linf"){
      uv_test_statistic <- apply(mv_test_statistic, c(2, 3), function(vec) max(abs(vec)))
      stagewise_p_values <- 1 - pchisq(uv_test_statistic^2, df = 1)^number_event_types
    } else {
      cat("Chosen norm of multivariate test statistic not available. Choose either \"l2\" or \"linf\"!")
    }

    # Set group sequential boundaries according to chosen method and number of analyses
    if(boundaries == "obf"){
      levels <- rpact::getDesignGroupSequential(kMax = number_of_analyses, alpha = alpha, sided = 2,
                                                typeOfDesign = "OF")$stageLevels * 2
    } else if(boundaries == "pocock"){
      levels <- rpact::getDesignGroupSequential(kMax = number_of_analyses, alpha = alpha, sided = 2,
                                                typeOfDesign = "P")$stageLevels * 2
    } else {
      cat("Chosen sequential boundary not available. Choose either \"obf\" or \"pocock\"!")
    }

    # Helper function for execution of inverse normal designs
    inverse_normal_combine <- function(p){
      k <- length(p)
      reweighting_factor <- sqrt(1/sum(weights[1:k]^2))
      current_weights <- weights[1:k] * reweighting_factor
      summed_z <- sum(current_weights*qnorm(1-p))
      return(1 - pnorm(summed_z))
    }

    # Compute p-values with inverse normal combination
    p_cum <- matrix(0, nrow = current_analysis, ncol = number_accrual_durations)

    for(i in 1:current_analysis){
      # Need to apply matrix function for the special case i=1, otherwise this object is not recognized as a matrix
      p_cum[i, ] <- apply(matrix(stagewise_p_values[1:i, ], nrow = i, ncol = number_accrual_durations), 2, inverse_normal_combine)
    }

    # Determine rejection
    rejection <- p_cum <= matrix(rep(levels[1:current_analysis], number_accrual_durations), ncol = number_accrual_durations)
    colnames(rejection) <- paste("a=", accrual_durations, sep = "")

    # Helper function for following use of 'apply'
    determine_overall_decision <- function(seq_decision){

      if(current_analysis == number_of_analyses & sum(seq_decision) == 0){
        decision <- "Accept H0"
        rejection_stage <- NA
      } else if(sum(seq_decision) >= 1){
        rejection_stage <- min(which(seq_decision == 1))
        decision <- paste("Reject H0 (in stage ", rejection_stage, ")", sep = "")
      } else {
        decision <- "Continue trial"
        rejection_stage <- NA
      }

      return(c(decision, rejection_stage))

    }

    decision_summary <- apply(rejection, 2, determine_overall_decision)

    decision <- decision_summary[1, ]
    rejection_stage <- decision_summary[2, ]

    output <- list(raw_martingale = mv_martingale,
                   stagewise_p_values = stagewise_p_values,
                   decision = decision,
                   rejection_stage = rejection_stage)

    return(output)
  }
}
