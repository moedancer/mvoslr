#' Helper function to create cumulative hazard function for exponential distribution
#'
#' @param rate Rate of exponential distribution
#'
#' @return conditional hazard function
#'
#' @keywords internal
#'
#' @examples
#' rate <- 1
#' reference_ch_fct <- mvoslr:::get_cum_haz_fct_exp(rate)
get_cum_haz_fct_exp <- function(rate){

  rate <- unname(rate)
  ch_fct <- function(t) rate*t
  return(ch_fct)

}

#' Helper function to create cumulative hazard function for exponential distribution
#'
#' @param parameters Vector of length two containing shape and scale parameter
#'
#' @return Conditional hazard function
#'
#' @keywords internal
#'
#' @examples
#' shape <- 1.5
#' scale <- 1
#' reference_ch_fct <- mvoslr:::get_cum_haz_fct_weibull(c(shape, scale))
get_cum_haz_fct_weibull <- function(parameters){

  parameters <- unname(parameters)
  ch_fct <- function(t) (t/parameters[2])^parameters[1]
  return(ch_fct)

}

#' Discretize list of functions for given time window and given number of steps
#'
#' @param function_list list of functions that shall be discretized
#' @param max_argument size of time window
#' @param time_steps number of steps
#'
#' @return Data frame containing values of functions at chosen input values
#'
#' @keywords internal
#'
#' @examples
#' f1 <- function(t) t^1.1
#' f2 <- function(t) t^0.9
#' f3 <- function(t) t^1.2
#' function_list_example <- list(f1, f2, f3)
#' mvoslr:::discretize_functions(function_list_example, max_argument = 10, time_steps = 1000)
discretize_functions <- function(function_list, max_argument, time_steps = 100){

  args <- seq(0, max_argument, max_argument/time_steps)

  values_frame <- data.frame(time = rep(args, length(function_list)),
                             Haz = unlist(lapply(function_list, function(x) x(args))),
                             trans = rep(1:length(function_list), each = time_steps + 1))

  return(values_frame)

}

#' Wrapper function for mstate::mssample with some adjustments
#'
#' @param transition_matrix Matrix of transitions between states as in mstate package
#' @param model_type Multi-state model is either Markov (\code{model_type = "M"}) or Semi-Markov (\code{model_type = "SM"})
#' @param cum_hazards_frame Data frame containing values of hazard functions of each transition
#' @param sample_size Number of subjects to simulate
#'
#' @return Data frame with results of simulation
#'
#' @keywords internal
#'
#' @examples
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' disc_cum_hazards_frame <- mvoslr:::discretize_functions(cum_hazards_example,
#'                                                         max_argument = 10, time_steps = 1000)
#' sample_size_example <- 100
#' mvoslr:::simulate_msm(transition_matrix = tmat_example, model_type = "M",
#'                       cum_hazards_frame = disc_cum_hazards_frame, sample_size = 10)
simulate_msm <- function(transition_matrix, model_type, cum_hazards_frame, sample_size){

  # Simulate multi-state data with mstate
  # Simulation might fail for unknown reasons, repeat simulation function until it does not fail
  successful_simulation <- FALSE
  while(!successful_simulation){
    if(model_type == "M"){
      sim_data <- try(mstate::mssample(Haz = cum_hazards_frame,
                                       trans = transition_matrix,
                                       M = sample_size,
                                       output = "data"), silent = TRUE)
    } else if(model_type == "SM"){
      sim_data <- try(mstate::mssample(Haz = cum_hazards_frame,
                                       trans = transition_matrix,
                                       M = sample_size,
                                       clock = "reset",
                                       output = "data"), silent = TRUE)
    }
    if (class(sim_data) != "try-error") successful_simulation <- TRUE
  }

  sim_frame <- msdata_to_df(sim_data)

  return(sim_frame)

}

#' Transform raw multi-state model data to trial data using simple censoring mechanism
#'
#' @param msm_data Data frame containing information about any (possible) transition with (at least) the following columns:
#' \itemize{
#'   \item id - ID of the patient to whom this entry belongs (requires continuous numbering starting with 1)
#'   \item Tstart - Start of the observation period for this transition
#'   \item Tstop - Stop of the observation period for this transition
#'   \item status - Transition could (1) or could not be (0) observed
#'   \item (recruitment_date) - If a recruitment date is supplied, it will be used to create corresponding censoring
#' }
#' @param accrual_duration Accrual duration of the trial to be simulated
#' @param follow_up_duration Duration of the follow-up period of the trial to be simulated
#'
#' @return Data frame with results of simulation
#'
#' @keywords internal
#'
#' @examples
#' raw_msm_data_example <- data.frame(id = c(1,1,2,2,2,3,3,4,4,4), Tstart = c(0,0,0,0,0.5,0,0,0,0,1),
#'                                    Tstop = c(0.8,0.8,0.5,0.5,1.2,1.3,1.3,1,1,1.1),
#'                                    status = c(0,1,1,0,1,0,0,1,0,0))
#' mvoslr:::msm_to_trial_data(raw_msm_data_example, accrual_duration = 1, follow_up_duration = 0.5)
msm_to_trial_data <- function(msm_data, accrual_duration, follow_up_duration){

  sample_size <- length(unique(msm_data$id))

  if(is.null(msm_data$recruitment_date)){
    recruitment_dates <- runif(n = sample_size,
                               min = 0,
                               max = accrual_duration)

    ids_occurences <- table(msm_data$id)
    msm_data$recruitment_date <- rep(recruitment_dates,
                                     ids_occurences)
  } else {
    msm_data <- msm_data[which(msm_data$recruitment_date <= accrual_duration),]
  }

  msm_data$censoring_date <- accrual_duration + follow_up_duration - msm_data$recruitment_date
  msm_data$status <- msm_data$status * (msm_data$Tstop <= msm_data$censoring_date)

  # Adapt end and duration of observation period according to censoring
  msm_data <- msm_data[which(msm_data$Tstart < msm_data$censoring_date),]
  msm_data$Tstop <- pmin(msm_data$Tstop, msm_data$censoring_date)
  msm_data$duration <- msm_data$Tstop - msm_data$Tstart

  msm_data$censoring_date <- NULL

  return(msm_data)

}

#' Transform data of type msdata from mstate package to a data frame
#'
#' @param msdata Data set of type msdata from \code{mstate} package
#'
#' @return data frame with column names corresponding to the columns names of msdata objects
#'
#' @keywords internal
#'
#' @examples
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' model_type_example <- "SM"
#' cumhaz_12_example <- function(t) t^1.1
#' cumhaz_13_example <- function(t) t^1.2
#' cumhaz_23_example <- function(t) t^0.9
#' cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
#' disc_cum_hazards_frame <- mvoslr:::discretize_functions(cum_hazards_example,
#'                                                         max_argument = 10, time_steps = 1000)
#' sim_data <- mstate::mssample(Haz = disc_cum_hazards_frame,
#'                              trans = tmat_example,
#'                              M = 10,
#'                              output = "data")
#' sim_frame <- msdata_to_df(sim_data)
msdata_to_df <- function(msdata){

  ms_names <- names(msdata)
  attributes(msdata) <- NULL

  ms_df <- data.frame(msdata)
  colnames(ms_df) <- ms_names

  return(ms_df)

}
