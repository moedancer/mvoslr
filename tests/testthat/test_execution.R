library(mvoslr)
library(mstate)

test_that("result of execution function of multivariate log-rank test coincides
          with calculation by hand for markovian illness-death model", {
  tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
  model_type_example <- "M"
  cumhaz_12_example <- function(t) t^2
  cumhaz_13_example <- function(t) t
  cumhaz_23_example <- function(t) sqrt(t)
  cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
  reference_model_example <- new_reference_model(transition_matrix = tmat_example,
                                                 intensities = cum_hazards_example,
                                                 type = model_type_example)
  analysis_dates_example <- c(3, 6)
  events_example <- list(c(2,3), c(3))
  names(events_example) <- c("PFS", "OS")
  msm_data <- data.frame(id               = c(1,1,1,2,  2,  3,3,3,4,4,5,  5,  5,  6,6,7,  7,  7,  8,8),
                         Tstart           = c(0,0,2,0,  0,  0,0,2,0,0,0,  0,  1.5,0,0,0,  0,  1.5,0,0),
                         Tstop            = c(2,2,4,1.5,1.5,2,2,3,3,3,1.5,1.5,3.5,2,2,1.5,1.5,2.5,2,2),
                         duration         = c(2,2,2,1.5,1.5,2,2,1,3,3,1.5,1.5,2,  2,2,1.5,1.5,1,  2,2),
                         from             = c(1,1,2,1,  1,  1,1,2,1,1,1,  1,  2,  1,1,1,  1,  2,  1,1),
                         to               = c(2,3,3,2,  3,  2,3,3,2,3,2,  3,  3,  2,3,2,  3,  3,  2,3),
                         status           = c(1,0,1,0,  1,  1,0,0,0,0,1,  0,  0,  0,1,1,  0,  0,  0,0),
                         trans            = c(1,2,3,1,  2,  1,2,3,1,2,1,  2,  3,  1,2,1,  2,  3,  1,2),
                         recruitment_date = c(0,0,0,0.5,0.5,1,1,1,2,2,2.5,2.5,2.5,3,3,3.5,3.5,3.5,4,4))

  # Results of calculation by hand
  N_1_PFS   <- 3
  N_1_OS    <- 1
  N_1_joint <- 1
  N_2_PFS   <- 3
  N_2_OS    <- 2
  N_2_joint <- 1
  A_1_PFS   <- 16+1.5^2+0.5^2
  A_1_OS    <- 7+3^0.5-2^0.5
  A_1_joint <- 7
  A_2_PFS   <- 24.5+2*1.5^2-0.5^2
  A_2_OS    <- 8.5+4^0.5-2^0.5+3.5^0.5+2.5^0.5-2*1.5^0.5
  A_2_joint <- 8.5

  cov_mat_1 <- 0.5 * matrix(c(N_1_PFS, N_1_joint, N_1_joint, N_1_OS), nrow = 2) +
    0.5 * matrix(c(A_1_PFS, A_1_joint, A_1_joint, A_1_OS), nrow = 2)
  cov_mat_2 <- 0.5 * matrix(c(N_2_PFS, N_2_joint, N_2_joint, N_2_OS), nrow = 2) +
    0.5 * matrix(c(A_2_PFS, A_2_joint, A_2_joint, A_2_OS), nrow = 2)

  std_mat_1 <- solve(t(chol(cov_mat_1)))
  std_mat_2 <- solve(t(chol(cov_mat_2)))

  # Result of R function
  result <- execution_mvoslr(msm_data = msm_data,
                             analysis_dates = analysis_dates_example,
                             current_analysis = 2,
                             reference_model = reference_model_example,
                             events = events_example)

  # Equality of of raw multivariate process
  expect_equal(unname(result$raw_martingale[,1]),
               c(N_1_PFS - A_1_PFS, N_1_OS - A_1_OS))
  expect_equal(unname(result$raw_martingale[,2]),
               c(N_1_PFS + N_2_PFS - (A_1_PFS + A_2_PFS),
                 N_1_OS + N_2_OS - (A_1_OS + A_2_OS)))

  # Equality of variance estimates
  expect_equal(unname(result$covariation_matrices[[1]]),
               cov_mat_1)
  expect_equal(unname(result$covariation_matrices[[2]]),
               cov_mat_1 + cov_mat_2)

  # Equality of standardised increments
  expect_equal(unname(result$multivariate_test_statistics[,1]),
               as.vector(std_mat_1 %*% c(N_1_PFS - A_1_PFS, N_1_OS - A_1_OS)))
  expect_equal(unname(result$multivariate_test_statistics[,2]),
               as.vector(std_mat_2 %*% c(N_2_PFS - A_2_PFS, N_2_OS - A_2_OS)))

  # Equality of univariate stagewise test statistics
  expect_equal(unname(result$stagewise_test_statistics[1]),
               sqrt(sum((std_mat_1 %*% c(N_1_PFS - A_1_PFS, N_1_OS - A_1_OS))^2)))
  expect_equal(unname(result$stagewise_test_statistics[2]),
               sqrt(sum((std_mat_2 %*% c(N_2_PFS - A_2_PFS, N_2_OS - A_2_OS))^2)))

  # Equality of stagewise p-values
  expect_equal(unname(result$stagewise_p_values[1]),
               1 - pchisq(sum((std_mat_1 %*% c(N_1_PFS - A_1_PFS, N_1_OS - A_1_OS))^2), df = 2))
  expect_equal(unname(result$stagewise_p_values[2]),
               1 - pchisq(sum((std_mat_2 %*% c(N_2_PFS - A_2_PFS, N_2_OS - A_2_OS))^2), df = 2))

})


test_that("result of execution function of multivariate log-rank test coincides
          with calculation by hand for markovian illness-death model", {
  tmat_example <- transMat(x = list(c(2,3,4),c(3,4),c(2,4),c()), names = c("h", "eff", "tox", "d"))
  model_type_example <- "SM"
  cumhaz_12_example <- function(t) t
  cumhaz_13_example <- function(t) t
  cumhaz_14_example <- function(t) t
  cumhaz_23_example <- function(t) t^2
  cumhaz_24_example <- function(t) t^2
  cumhaz_32_example <- function(t) t^3
  cumhaz_34_example <- function(t) t^3
  cum_hazards_example <- list(cumhaz_12_example,
                              cumhaz_13_example,
                              cumhaz_14_example,
                              cumhaz_23_example,
                              cumhaz_24_example,
                              cumhaz_32_example,
                              cumhaz_34_example)
  reference_model_example <- new_reference_model(transition_matrix = tmat_example,
                                                 intensities = cum_hazards_example,
                                                 type = model_type_example)
  analysis_dates_example <- c(2,4,6)
  events_example <- list(c(2,4), c(3,4))
  names(events_example) <- c("Eff", "Tox")
  msm_data <- data.frame(id               = c(1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,5,  5,  5,  5,  5,  6,6,6,6,6,7,7,7,8,  8,  8,  8,  8  ),
                         Tstart           = c(0,0,0,1,1,0,0,0,2,2,3,3,0,0,0,3,3,1,1,0,0,0,0,  0,  0,  0.5,0.5,0,0,0,1,1,0,0,0,0,  0,  0,  0.5,0.5),
                         Tstop            = c(1,1,1,3,3,2,2,2,3,3,5,5,1,1,1,5,5,3,3,2,2,2,0.5,0.5,0.5,1,  1,  1,1,1,2,2,2,2,2,0.5,0.5,0.5,1,  1  ),
                         duration         = c(1,1,1,2,2,2,2,2,1,1,2,2,1,1,1,2,2,2,2,2,2,2,0.5,0.5,0.5,0.5,0.5,1,1,1,1,1,2,2,2,0.5,0.5,0.5,0.5,0.5),
                         from             = c(1,1,1,2,2,1,1,1,2,2,3,3,1,1,1,2,2,3,3,1,1,1,1,  1,  1,  2,  2,  1,1,1,3,3,1,1,1,1,  1,  1,  3,  3  ),
                         to               = c(2,3,4,3,4,2,3,4,3,4,2,4,2,3,4,3,4,2,4,2,3,4,2,  3,  4,  3,  4,  2,3,4,2,4,2,3,4,2,  3,  4,  2,  4  ),
                         status           = c(1,0,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,1,0,0,0,1,1,  0,  0,  0,  1,  0,1,0,0,1,0,0,0,0,  1,  0,  0,  0  ),
                         trans            = c(1,2,3,4,5,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,1,  2,  3,  4,  5,  1,2,3,6,7,1,2,3,1,  2,  3,  6,  7  ),
                         recruitment_date = c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,3,  3,  3,  3,  3,  3,3,3,3,3,4,4,4,5,  5,  5,  5,  5  ))

  # Results of calculation by hand
  N_1_Eff   <- 2
  N_1_Tox   <- 1
  N_1_joint <- 0
  N_2_Eff   <- 3
  N_2_Tox   <- 5
  N_2_joint <- 1
  N_3_Eff   <- 1
  N_3_Tox   <- 1
  N_3_joint <- 0
  A_1_Eff   <- 8
  A_1_Tox   <- 10
  A_1_joint <- 4
  A_2_Eff   <- 23
  A_2_Tox   <- 15.5
  A_2_joint <- 3.5
  A_3_Eff   <- 7.25
  A_3_Tox   <- 5
  A_3_joint <- 2.5


  cov_mat_1 <- 0.5 * matrix(c(N_1_Eff, N_1_joint, N_1_joint, N_1_Tox), nrow = 2) +
    0.5 * matrix(c(A_1_Eff, A_1_joint, A_1_joint, A_1_Tox), nrow = 2)
  cov_mat_2 <- 0.5 * matrix(c(N_2_Eff, N_2_joint, N_2_joint, N_2_Tox), nrow = 2) +
    0.5 * matrix(c(A_2_Eff, A_2_joint, A_2_joint, A_2_Tox), nrow = 2)
  cov_mat_3 <- 0.5 * matrix(c(N_3_Eff, N_3_joint, N_3_joint, N_3_Tox), nrow = 2) +
    0.5 * matrix(c(A_3_Eff, A_3_joint, A_3_joint, A_3_Tox), nrow = 2)

  std_mat_1 <- solve(t(chol(cov_mat_1)))
  std_mat_2 <- solve(t(chol(cov_mat_2)))
  std_mat_3 <- solve(t(chol(cov_mat_3)))

  # Result of R function
  result <- execution_mvoslr(msm_data = msm_data,
                             analysis_dates = analysis_dates_example,
                             current_analysis = 3,
                             reference_model = reference_model_example,
                             events = events_example)

  # Equality of of raw multivariate process
  expect_equal(unname(result$raw_martingale[,1]),
               c(N_1_Eff - A_1_Eff, N_1_Tox - A_1_Tox))
  expect_equal(unname(result$raw_martingale[,2]),
               c(N_1_Eff + N_2_Eff - (A_1_Eff + A_2_Eff),
                 N_1_Tox + N_2_Tox - (A_1_Tox + A_2_Tox)))
  expect_equal(unname(result$raw_martingale[,3]),
               c(N_1_Eff + N_2_Eff + N_3_Eff - (A_1_Eff + A_2_Eff + A_3_Eff),
                 N_1_Tox + N_2_Tox + N_3_Tox - (A_1_Tox + A_2_Tox + A_3_Tox)))

  # Equality of variance estimates
  expect_equal(unname(result$covariation_matrices[[1]]),
               cov_mat_1)
  expect_equal(unname(result$covariation_matrices[[2]]),
               cov_mat_1 + cov_mat_2)
  expect_equal(unname(result$covariation_matrices[[3]]),
               cov_mat_1 + cov_mat_2 + cov_mat_3)

  # Equality of standardised increments
  expect_equal(unname(result$multivariate_test_statistics[,1]),
               as.vector(std_mat_1 %*% c(N_1_Eff - A_1_Eff, N_1_Tox - A_1_Tox)))
  expect_equal(unname(result$multivariate_test_statistics[,2]),
               as.vector(std_mat_2 %*% c(N_2_Eff - A_2_Eff, N_2_Tox - A_2_Tox)))
  expect_equal(unname(result$multivariate_test_statistics[,3]),
               as.vector(std_mat_3 %*% c(N_3_Eff - A_3_Eff, N_3_Tox - A_3_Tox)))

  # Equality of univariate stagewise test statistics
  expect_equal(unname(result$stagewise_test_statistics[1]),
               sqrt(sum((std_mat_1 %*% c(N_1_Eff - A_1_Eff, N_1_Tox - A_1_Tox))^2)))
  expect_equal(unname(result$stagewise_test_statistics[2]),
               sqrt(sum((std_mat_2 %*% c(N_2_Eff - A_2_Eff, N_2_Tox - A_2_Tox))^2)))
  expect_equal(unname(result$stagewise_test_statistics[3]),
               sqrt(sum((std_mat_3 %*% c(N_3_Eff - A_3_Eff, N_3_Tox - A_3_Tox))^2)))

  # Equality of stagewise p-values
  expect_equal(unname(result$stagewise_p_values[1]),
               1 - pchisq(sum((std_mat_1 %*% c(N_1_Eff - A_1_Eff, N_1_Tox - A_1_Tox))^2), df = 2))
  expect_equal(unname(result$stagewise_p_values[2]),
               1 - pchisq(sum((std_mat_2 %*% c(N_2_Eff - A_2_Eff, N_2_Tox - A_2_Tox))^2), df = 2))
  expect_equal(unname(result$stagewise_p_values[3]),
               1 - pchisq(sum((std_mat_3 %*% c(N_3_Eff - A_3_Eff, N_3_Tox - A_3_Tox))^2), df = 2))

})
