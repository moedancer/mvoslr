#' Execution of interim and final analyses of multivariate group-sequential one-sample log-rank test
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
#' @param analysis_dates Vector of calendar dates of analyses
#' @param current_analysis Number of the current analysis
#' @param reference_model Specification of the reference model against which the new data is tested. Should be an object of class "reference_model"
#' @param events List of (composite) events that shall be investigated
#' @param norm Use either \eqn{L^2}-norm (\code{norm = "l2"}, default value) or \eqn{L^\infty}-norm (\code{norm = "linf"}) of vector of test statistics to compute stagewise p-values
#' @param boundaries Use either O'Brien-Fleming'S (\code{boundaries = "obf"}, default value) or Pocock's (\code{boundaries = "pocock"}) sequential decision boundaries
#' @param alpha Choose type I error rate (default value = 0.05)
#' @param weights Choose weights for inverse normal combination of stagewise p-values. Sum of squared values needs to sum up to 1.
#'
#' @return Object of class "mvoslr_result"
#'
#' @import stats
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
#' reference_model_example <- new_reference_model(transition_matrix = tmat_example,
#'                                                intensities = cum_hazards_example,
#'                                                type = model_type_example)
#' analysis_dates_example <- c(1, 2)
#' events_example <- list(c(2,3), c(3))
#' names(events_example) <- c("PFS", "OS")
#' #Fictional data from observations in multi-state model with one dead patient without illness
#' #  (\code{id =1}), one dead patient with illness (\code{id =2}), one healthy, censored patient
#' #  (\code{id =3}) and one ill, censored patient (\code{id =4})
#' msm_data_example <- data.frame(id = c(1,1,2,2,2,3,3,4,4,4), Tstart = c(0,0,0,0,0.5,0,0,0,0,1),
#'                                Tstop = c(0.8,0.8,0.5,0.5,1.2,1.3,1.3,1,1,1.1),
#'                                duration = c(0.8,0.8,0.5,0.5,0.7,1.3,1.3,1,1,0.1),
#'                                from = c(1,1,1,1,2,1,1,1,1,2), to = c(2,3,2,3,3,2,3,2,3,3),
#'                                status = c(0,1,1,0,1,0,0,1,0,0), trans = c(1,2,1,2,3,1,2,1,2,3),
#'                                recruitment_date = c(0,0,0.3,0.3,0.3,0.7,0.7,0.9,0.9,0.9))
#' execution_mvoslr(msm_data_example, analysis_dates = analysis_dates_example, current_analysis = 2,
#'                  reference_model = reference_model_example, events = events_example)
execution_mvoslr <- function(msm_data, analysis_dates, current_analysis = NULL, reference_model,
                             events, norm = "l2", boundaries = "obf", alpha = 0.05, weights = NULL){

  # Unpack information from reference model
  transition_matrix <- reference_model$transition_matrix
  cum_hazard_functions <- reference_model$intensities
  model_type <- attributes(reference_model)$type

  if(!model_type %in% c("M", "SM")){
    cat("Parameter model needs to be one of \"M\" (Markov model) or \"SM\" (Semi-Markov model)!" )
  } else if(!norm %in% c("l2", "linf")) {
    cat("Chosen norm need to be either \"l2\" or \"linf\"!")
  } else if(!boundaries %in% c("obf", "pocock")) {
    cat("Sequential boundaries needs to be one of \"obf\" or \"pocock\"!")
  } else if(!is.null(weights) & abs(1 - sum(weights^2)) > 0.00001){
    cat("Squared weights do not sum up to one!")
  } else if(!is.null(weights) & length(weights)!=length(analysis_dates)){
    cat("Length of weight vector does not correspond to number of analyses!")
  } else {

    number_event_types <- length(events)
    number_of_analyses <- length(analysis_dates)

    transitions <- mstate::to.trans2(transition_matrix)
    number_of_trans <- dim(transitions)[1]

    if(is.null(weights)) weights <- rep(1/sqrt(number_of_analyses), number_of_analyses)

    # Current analysis is set to overall number of analyses if not specified otherwise
    if(is.null(current_analysis)) current_analysis <- length(analysis_dates)

    #create arrays to store cumulative and incremental events and accumulated hazards
    events_count  <- array(0, c(number_event_types, number_event_types, current_analysis))
    events_hazard <- array(0, c(number_event_types, number_event_types, current_analysis))

    events_count_inc  <- array(0, c(number_event_types, number_event_types, current_analysis))
    events_hazard_inc <- array(0, c(number_event_types, number_event_types, current_analysis))

    #create matrices/vectors for martingale, test statistics and p-values
    mv_martingale <- matrix(0, nrow = number_event_types, ncol = current_analysis)
    mv_martingale_inc <- matrix(0, nrow = number_event_types, ncol = current_analysis)
    mv_test_statistic <- matrix(0, nrow = number_event_types, ncol = current_analysis)
    uv_test_statistic <- rep(0, current_analysis)

    if(!is.null(names(events))){
      rownames(mv_martingale) <- names(events)
      rownames(mv_test_statistic) <- names(events)
      dimnames(events_count)[[1]] <- dimnames(events_count)[[2]] <- names(events)
      dimnames(events_hazard)[[1]] <- dimnames(events_hazard)[[2]] <- names(events)
      dimnames(events_count_inc)[[1]] <- dimnames(events_count_inc)[[2]] <- names(events)
      dimnames(events_hazard_inc)[[1]] <- dimnames(events_hazard_inc)[[2]] <- names(events)
    }

    for(number_of_analysis in 1:current_analysis){

      analysis_date <- analysis_dates[number_of_analysis]

      msm_data_temp <- msm_data

      msm_data_temp$censoring_date <- analysis_date - msm_data$recruitment_date

      #exclude patients which were not recruited at date of analysis
      msm_data_temp <- msm_data_temp[which(msm_data_temp$censoring_date >= 0), ]

      #exclude observations starting after censoring
      msm_data_temp <- msm_data_temp[which(msm_data_temp$Tstart <= msm_data_temp$censoring_date), ]

      #adapt status of transitions according to censoring
      msm_data_temp$status <- msm_data_temp$status * (msm_data_temp$Tstop <= msm_data_temp$censoring_date)

      #adapt end and duration of observation period according to censoring
      msm_data_temp$Tstop <- pmin(msm_data_temp$Tstop, msm_data_temp$censoring_date)
      msm_data_temp$duration <- msm_data_temp$Tstop - msm_data_temp$Tstart

      ### Step 8: Calculate accumulated hazards for each transition

      #create column to enter accumulated hazards
      msm_data_temp$acc_haz <- rep(0, dim(msm_data_temp)[1])

      for(i in 1:number_of_trans){

        current_cumhaz <- cum_hazard_functions[[i]]
        current_start <- transitions$from[i]
        current_end <- transitions$to[i]

        current_transitions <- msm_data_temp[which(msm_data_temp$from == current_start &
                                                     msm_data_temp$to == current_end), ]

        if(model_type == "M"){
          current_acc_hazards <- current_cumhaz(current_transitions$Tstop) -
            current_cumhaz(current_transitions$Tstart)
        } else if(model_type == "SM"){
          current_acc_hazards <- current_cumhaz(current_transitions$duration)
        }

        msm_data_temp$acc_haz[which(msm_data_temp$from == current_start &
                                      msm_data_temp$to == current_end)] <- current_acc_hazards

      }

      for(i in 1:number_event_types){

        for(j in i:number_event_types){

          event_1 <- events[[i]]
          event_2 <- events[[j]]
          union_event <- sort(unique(c(event_1, event_2)))
          intersect_event <- intersect(event_1, event_2)

          if(length(intersect_event) > 0){

            #extract relevant information for current event
            relevant_obs <- msm_data_temp[which(msm_data_temp$to %in% union_event), ]
            true_events <- relevant_obs[which(relevant_obs$status == 1), ]

            #compute event time (add dummy rows if event did not take place)
            event_times <- aggregate(Tstop ~ id, true_events, function(x) min(x))
            all_patients <- unique(msm_data_temp$id)

            eventless_patients <- all_patients[which(!(all_patients %in% event_times$id))]
            number_eventless_patients <- length(eventless_patients)

            if(number_eventless_patients > 0){
              dummy_rows <- data.frame(id = eventless_patients, Tstop = rep(Inf, number_eventless_patients))
              event_times <- rbind(event_times, dummy_rows)
            }

            event_times <- event_times[order(event_times$id),]

            if(!all(event_times$id == all_patients)) print("HELP!")

            #extract observations until occurence of event
            event_times_rep <-  event_times$Tstop[relevant_obs$id]
            relevant_obs <- relevant_obs[which(relevant_obs$Tstop <= event_times_rep &
                                                 relevant_obs$to %in% intersect_event), ]

            sum_intersect_events <- sum(relevant_obs$status)
            sum_accumulated_hazard <- sum(relevant_obs$acc_haz)

            events_count[i, j, number_of_analysis]  <- events_count[j, i, number_of_analysis]  <- sum_intersect_events
            events_hazard[i, j, number_of_analysis] <- events_hazard[j, i, number_of_analysis] <- sum_accumulated_hazard

          }
        }
      }

      mv_martingale[,number_of_analysis] <- diag(events_count[,,number_of_analysis]) - diag(events_hazard[,,number_of_analysis])

      #compute increments of event count and cumulative hazard data
      if(number_of_analysis == 1){

        events_count_inc[,,1] <- events_count[,,1]
        events_hazard_inc[,,1] <- events_hazard[,,1]
        mv_martingale_inc[,1] <- mv_martingale[,1]

      } else {

        events_count_inc[,,number_of_analysis] <- events_count[,,number_of_analysis] - events_count[,,number_of_analysis-1]
        events_hazard_inc[,,number_of_analysis] <- events_hazard[,,number_of_analysis] - events_hazard[,,number_of_analysis-1]
        mv_martingale_inc[,number_of_analysis] <- mv_martingale[,number_of_analysis] - mv_martingale[,number_of_analysis-1]

      }

    }

    covariance_estimate <- 0.5 * events_count + 0.5 * events_hazard
    covariance_estimate_inc <- 0.5 * events_count_inc + 0.5 * events_hazard_inc

    standardisation_matrices_raw <- apply(covariance_estimate_inc, 3, function(matrix) solve(t(chol(matrix))))
    standardisation_matrices <- array(standardisation_matrices_raw, dim = c(number_event_types, number_event_types, current_analysis))

    mv_test_statistic <- matrix(0, nrow = number_event_types, ncol = current_analysis)
    for(i in 1:current_analysis) mv_test_statistic[,i] <- standardisation_matrices[,,i] %*% mv_martingale_inc[,i]

    if(!is.null(names(events))){
      rownames(mv_test_statistic) <- names(events)
    }

    if(norm == "l2"){
      uv_test_statistic <- apply(mv_test_statistic, 2, function(vec) sqrt(sum(vec^2)))
      stagewise_p_values <- 1 - pchisq(uv_test_statistic^2, df = number_event_types)
    } else if(norm == "linf"){
      uv_test_statistic <- apply(mv_test_statistic, 2, function(vec) max(abs(vec)))
      stagewise_p_values <- 1 - pchisq(uv_test_statistic^2, df = 1)^number_event_types
    } else {
      cat("Chosen norm of multivariate test statistic not available. Choose either \"l2\" or \"linf\"!")
    }

    if(boundaries == "obf"){
      levels <- rpact::getDesignGroupSequential(kMax = number_of_analyses, alpha = alpha, sided = 2,
                                                typeOfDesign = "OF")$stageLevels * 2
    } else if(boundaries == "pocock"){
      levels <- rpact::getDesignGroupSequential(kMax = number_of_analyses, alpha = alpha, sided = 2,
                                                typeOfDesign = "P")$stageLevels * 2
    } else {
      cat("Chosen sequential boundary not available. Choose either \"obf\" or \"pocock\"!")
    }

    # Define 'local' version of inverse normal combination function with weights specified above
    inverse_normal_combine_loc <- function(p){
      inverse_normal_combine(p, weights)
    }

    p_cum <- rep(0, current_analysis)

    # Combine p_values for each stage and determine decision
    for(i in 1:current_analysis){
      p_cum[i] <- inverse_normal_combine_loc(stagewise_p_values[1:i])
    }

    rejection <- p_cum <= levels[1:current_analysis]

    # Rename covariance estimate if possible
    if(!is.null(names(events))){
      dimnames(covariance_estimate)[[1]] <- dimnames(covariance_estimate)[[2]]  <- names(events)
    }

    # Extract rejection stage or set it to NA if no rejection occure
    if(length(which(rejection)) == 0) {rejection_stage <- NA} else {rejection_stage <- min(which(rejection))}

    result <- new_mvoslr_result(raw_process = mv_martingale, covariance_matrices = covariance_estimate,
                                multivariate_test_statistics = mv_test_statistic,
                                univariate_test_statistics = uv_test_statistic,
                                stagewise_p_values = stagewise_p_values, rejection_stage = rejection_stage,
                                remaining_analyses = as.integer(length(analysis_dates) - current_analysis),
                                vector_norm = norm)

    return(result)
  }
}
