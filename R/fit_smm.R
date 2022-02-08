#' Fitting of Semi-Markov model with Weibull transition intensities to (historic) data from multi-state model
#'
#' @param msm_data data frame containing information about any (possible) transition of patients in the study with the following columns:
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
#' @param transition_matrix matrix of transitions between states as in mstate package
#' @param joint_shape boolean variable indicating whether shape parameter is the same for all transitions (as in Li and Zhang (2015)) or not
#'
#' @return List of 2:
#' \itemize{
#'   \item parameters - Estimated parameters for each transition
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
#' fit_smm(msm_data_example, transition_matrix = tmat_example)
#'
fit_smm <- function(msm_data, transition_matrix, joint_shape = TRUE){

  transitions <- mstate::to.trans2(transition_matrix)
  if(!all(sort(unique(msm_data$trans)) == transitions$transno)){
    cat("Not all transition intensities can be estimated. No data for transition(s) ",
        setdiff(transitions$transno, unique(msm_data$trans)), ".", sep = "")
  }

  # Prepare matrix to store parameters for each transition
  parameter_estimates <- matrix(NA, nrow = 2, ncol = length(transitions$transno))
  rownames(parameter_estimates) <- c("shape", "scale")
  colnames(parameter_estimates) <- transitions$transno

  if(!joint_shape){

      for(trans in transitions$transno){

      i_trans <- which(transitions$transno == trans)

      # Prepare data for fitdistcens function
      temp_data <- data.frame(left = msm_data$duration[which(msm_data$trans == trans)],
                              right = ifelse(msm_data$status[which(msm_data$trans == trans)] == 1,
                                             msm_data$duration[which(msm_data$trans == trans)],
                                             NA))
      est <- fitdistrplus::fitdistcens(temp_data,
                                       distr = "weibull")

      parameter_estimates[, i_trans] <- est$estimate

    }

  } else {

    log_likelihood <- function(parameters){
      return(-sum(apply(X = msm_data, MARGIN = 1, FUN = log_lik_single_obs_smm,
                        shape = parameters[1],
                        scales = parameters[-1])))
    }

    optim_result <- optim(rep(1, 1 + length(transitions$transno)),
                          fn = log_likelihood)

    parameter_estimates["shape", ] <- rep(optim_result[[1]][1], length(transitions$transno))
    parameter_estimates["scale", ] <- optim_result[[1]][-1]

  }

  # Convert parameter collection to list for creation of cumulative hazard functions
  parameter_estimates_list <- lapply(1:ncol(parameter_estimates),
                                     function(col_number) parameter_estimates[, col_number])

  cum_hazard_functions <- lapply(parameter_estimates_list,
                                 FUN = get_cum_haz_fct_weibull)

  output <- list(parameter_estimates,
                 cum_hazard_functions)

  return(output)

}

#' Helper function to set up log likelihood function to estimate parameters in Semi-Markov model
#'
#' @param single_obs_data single row of data from multi-state model, column names
#' @param shape joint shape parameter of transition intensity functions
#' @param scales vector of scale parameters of transition intensity functions
#'
#' @return log likelihood for single observation
#'
#' @keywords internal
#'
#' @examples
#' msm_data_example <- data.frame(id = c(1,1,2,2,2,3,3,4,4,4), Tstart = c(0,0,0,0,0.5,0,0,0,0,1),
#'                                Tstop = c(0.8,0.8,0.5,0.5,1.2,1.3,1.3,1,1,1.1),
#'                                duration = c(0.8,0.8,0.5,0.5,0.7,1.3,1.3,1,1,0.1),
#'                                from = c(1,1,1,1,2,1,1,1,1,2), to = c(2,3,2,3,3,2,3,2,3,3),
#'                                status = c(0,1,1,0,1,0,0,1,0,0), trans = c(1,2,1,2,3,1,2,1,2,3),
#'                                recruitment_date = c(0,0,0.3,0.3,0.3,0.7,0.7,0.9,0.9,0.9))
#' mvoslr:::log_lik_single_obs_smm(msm_data_example[1,], shape = 1.5, scales = c(0.4, 0.2, 0.5))
log_lik_single_obs_smm <- function(single_obs_data, shape, scales){

  trans <- single_obs_data["trans"][[1]]
  status  <- single_obs_data["status"][[1]]
  duration <- single_obs_data["duration"][[1]]

  if(status == 1){
    log_lik <- log(dweibull(duration, shape = shape, scale = scales[trans]))
  } else {
    log_lik <- log(1 - pweibull(duration, shape = shape, scale = scales[trans]))
  }

  return(log_lik)

}
