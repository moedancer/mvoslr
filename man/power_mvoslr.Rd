% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_mvoslr.R
\name{power_mvoslr}
\alias{power_mvoslr}
\title{Estimation of power of group-sequential multivariate one-sample log-rank test via simulation}
\usage{
power_mvoslr(
  reference_model,
  events,
  analysis_dates,
  accrual_duration,
  sample_size,
  hazard_ratios = NULL,
  cum_hazard_functions_alternative = NULL,
  norm = "l2",
  boundaries = "obf",
  alpha = 0.05,
  weights = NULL,
  time_steps = 100,
  simulation_runs = 1000
)
}
\arguments{
\item{reference_model}{Specification of the reference model against which the new data is tested. Should be an object of class "reference_model".}

\item{events}{List of (composite) events that shall be investigated}

\item{analysis_dates}{Vector of calendar dates of analyses}

\item{accrual_duration}{Duration of accrual period}

\item{sample_size}{Sample size for which power shall be estimated}

\item{hazard_ratios}{Specification of alternative hypothesis. One can either choose a single hazard ratio which then applies
to each transition or a single hazard ratio for each transition in the model.}

\item{cum_hazard_functions_alternative}{If non-proportional transition intensities are anticipated, the transition intensities
for all transitions, need to be specified here. If this argument is specified, it
overrules the argument \code{hazard_ratios}.}

\item{norm}{Use either \eqn{L^2}-norm (\code{norm = "l2"}, default value) or \eqn{L^\infty}-norm (\code{norm = "linf"}) of vector of test statistics to compute stagewise p-values}

\item{boundaries}{Use either O'Brien-Fleming'S (\code{boundaries = "obf"}, default value) or Pocock's (\code{boundaries = "pocock"}) sequential decision boundaries}

\item{alpha}{Choose type I error rate (default value = 0.05)}

\item{weights}{Choose weights for inverse normal combination of stagewise p-values. Sum of squared values needs to sum up to 1.}

\item{time_steps}{As the multi-state model data will be simulated with the \code{mstate}-package, the transition intensities need
to be discretized. The number of time steps in which the intensities will be discretized can be specified here.
The time horizon over which the discretization happens is the calndar date of the last analysis.}

\item{simulation_runs}{The number of simulation runs to estimate the power can be specified here.}
}
\value{
Object of class "mvoslr_power_object"
}
\description{
Estimation of power of group-sequential multivariate one-sample log-rank test via simulation
}
\examples{
#Setup of reference multi-state model (here: simple illness-death model)
library(mstate)
tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
number_of_trans_example <- dim(to.trans2(tmat_example))[1]
model_type_example <- "SM"
cumhaz_12_example <- function(t) t^1.1
cumhaz_13_example <- function(t) t^1.2
cumhaz_23_example <- function(t) t^0.9
cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
reference_model_example <- new_reference_model(transition_matrix = tmat_example,
                                               intensities = cum_hazards_example,
                                               type = model_type_example)
analysis_dates_example <- c(1, 2)
events_example <- list(c(2,3), c(3))
names(events_example) <- c("PFS", "OS")
accrual_duration_example <- 1
sample_size_example <- 100
#In this example, the alternative is specified via separate hazard ratios for each transition
hazard_ratios_example <- c(1.4, 1.2, 1.35)
power_mvoslr(reference_model = reference_model_example, events = events_example,
             analysis_dates = analysis_dates_example,
             accrual_duration = accrual_duration_example, sample_size = sample_size_example,
             hazard_ratios = hazard_ratios_example, simulation_runs = 10)
}
