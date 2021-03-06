% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvoslr_result.R
\name{new_mvoslr_result}
\alias{new_mvoslr_result}
\title{Constructor for class "mvoslr_result" which is used to save and display results of multivariate one-sample log-rank tests}
\usage{
new_mvoslr_result(
  raw_process = matrix(),
  covariance_matrices = array(),
  multivariate_test_statistics = matrix(),
  univariate_test_statistics = numeric(),
  stagewise_p_values = numeric(),
  rejection_stage = integer(),
  remaining_analyses = integer(),
  vector_norm = character()
)
}
\arguments{
\item{raw_process}{Values of multivariate process at analysis dates}

\item{covariance_matrices}{Estimated covariance matrix for multivariate process at analysis date}

\item{multivariate_test_statistics}{Stagewise multivariate, standardised test-statistics}

\item{stagewise_p_values}{Stagewise p-values}

\item{rejection_stage}{Number of stage at which the null hypothesis was rejected (is \code{NA} if null hypothesis could not be rejected)}

\item{univariate_test_statistic}{Stagewise univariate test-statistics}
}
\value{
Object of class "mvoslr_result"
}
\description{
Constructor for class "mvoslr_result" which is used to save and display results of multivariate one-sample log-rank tests
}
\keyword{internal}
