#' Fitting of time-homogeneous Markov model to (historic) data from multi-state model
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
#' }
#' @param transition_matrix Matrix of transitions between states as in mstate package
#'
#' @return List of 2:
#' \itemize{
#'   \item parameter_estimates - Estimated parameters for each transition
#'   \item cum_hazard_function - List of cumulative hazard function to be passed to power analysis or execution of multivariate one-sample log-rank test
#' }
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' #Setup of reference multi-state model (here: simple illness-death model)
#' library(mstate)
#' tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
#' #Fictional data from observations in multi-state model with one dead patient without illness
#' #  (\code{id =1}), one dead patient with illness (\code{id =2}), one healthy, censored patient
#' #  (\code{id =3}) and one ill, censored patient (\code{id =4})
#' msm_data_example <- data.frame(id = c(1,1,2,2,2,3,3,4,4,4), Tstart = c(0,0,0,0,0.5,0,0,0,0,1),
#'                                Tstop = c(0.8,0.8,0.5,0.5,1.2,1.3,1.3,1,1,1.1),
#'                                duration = c(0.8,0.8,0.5,0.5,0.7,1.3,1.3,1,1,0.1),
#'                                from = c(1,1,1,1,2,1,1,1,1,2), to = c(2,3,2,3,3,2,3,2,3,3),
#'                                status = c(0,1,1,0,1,0,0,1,0,0), trans = c(1,2,1,2,3,1,2,1,2,3),
#'                                recruitment_date = c(0,0,0.3,0.3,0.3,0.7,0.7,0.9,0.9,0.9))
#' fit_thmm(msm_data_example, transition_matrix = tmat_example)
fit_thmm <- function(msm_data, transition_matrix){

  transitions <- mstate::to.trans2(transition_matrix)
  if(!all(sort(unique(msm_data$trans)) == transitions$transno)){
    cat("Not all transition intensities can be estimated. No data for transition(s) ",
        setdiff(transitions$transno, unique(msm_data$trans)), ".", sep = "")
  }

  num_trans <- dim(transitions)[1]

  # Aggregate time spent waiting for each transition and number of transitions
  agg_time   <- aggregate(duration ~ trans, data = msm_data, FUN = sum)
  agg_events <- aggregate(status ~ trans, data = msm_data, FUN = sum)

  agg_data <- merge(agg_time, agg_events, by = "trans")
  agg_data <- agg_data[order(agg_data$trans), ]

  # Compute estimated transition rates
  agg_data$rates <- agg_data$status/agg_data$duration

  # Extract required parameter estimates
  parameter_estimates <- agg_data[which(agg_data$trans %in% transitions$transno), "rates"]
  names(parameter_estimates) <- agg_data$trans[agg_data$trans %in% transitions$transno]

  # Transform parameters estimates to Weibull form
  parameter_estimates_weibull <- rbind(rep(1, num_trans), 1/parameter_estimates)
  rownames(parameter_estimates_weibull) <- c("shape", "scale")

  output_model <- reference_model_weibull(transition_matrix = transition_matrix,
                                          type = "M",
                                          parameters = parameter_estimates_weibull)

  return(output_model)

}
