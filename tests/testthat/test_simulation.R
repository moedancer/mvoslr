library(mstate)
library(mvoslr)

test_that("discretization function operates as desired", {
  my_function_1 <- function(t) t^2 + 1.5 * t * log(t + 1)
  my_function_2 <- function(t) sqrt(t) - t/4
  my_functions <- list(my_function_1, my_function_2)
  upper_bound <- 100
  steps <- 100

  result <- discretize_functions(my_functions, upper_bound, time_steps = steps)

  expect_equal(result$time[1], 0)
  expect_equal(result$time[101], 100)
  expect_equal(result$time[102], 0)

  expect_equal(result$trans[1], 1)
  expect_equal(result$trans[101], 1)
  expect_equal(result$trans[102], 2)

  expect_equal(result$Haz[1], 0)
  expect_equal(result$Haz[101], 100^2 + 1.5 * 100 * log(101))
  expect_equal(result$Haz[102], 0)
})

test_that("data is properly transformed into a data frame", {
  tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
  cumhaz_12_example <- function(t) t^1.2
  cumhaz_13_example <- function(t) t^1.3
  cumhaz_23_example <- function(t) t^0.4
  cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
  disc_cum_hazards_frame <- mvoslr:::discretize_functions(cum_hazards_example,
                                                          max_argument = 10, time_steps = 1000)
  sim_data <- mstate::mssample(Haz = disc_cum_hazards_frame,
                               trans = tmat_example,
                               M = 10,
                               output = "data")
  sim_frame <- msdata_to_df(sim_data)

  expect_equal(sim_frame$id, sim_data$id)
  expect_equal(sim_frame$Tstart, sim_data$Tstart)
})

test_that("output of simulation function is reproducible", {
  my_seed <- 100

  tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
  cumhaz_12_example <- function(t) t
  cumhaz_13_example <- function(t) t
  cumhaz_23_example <- function(t) t
  cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
  disc_cum_hazards_frame <- discretize_functions(cum_hazards_example,
                                                 max_argument = 10, time_steps = 1000)
  sample_size_example <- 100

  set.seed(my_seed)
  sim_data1 <- mssample(Haz = disc_cum_hazards_frame,
                        trans = tmat_example,
                        M = sample_size_example,
                        output = "data")
  sim_frame1 <- msdata_to_df(sim_data1)

  set.seed(my_seed)
  sim_frame2 <- simulate_msm(transition_matrix = tmat_example, model_type = "M",
                             cum_hazards_frame = disc_cum_hazards_frame, sample_size = sample_size_example)

  expect_equal(sim_frame1, sim_frame2)
})

test_that("simulated trial data has required properties", {
  tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
  model_type_example <- "SM"
  cumhaz_12_example <- function(t) t^1
  cumhaz_13_example <- function(t) t^2
  cumhaz_23_example <- function(t) t^0.5
  cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
  disc_cum_hazards_frame <- discretize_functions(cum_hazards_example,
                                                 max_argument = 10, time_steps = 1000)
  raw_frame <- simulate_msm(transition_matrix = tmat_example, model_type = "M",
                            cum_hazards_frame = disc_cum_hazards_frame, sample_size = 1000)

  a <- 4
  f <- 6

  trial_frame <- msm_to_trial_data(raw_frame, accrual_duration = a, follow_up_duration = f)

  expect_gte(min(trial_frame$recruitment_date), 0)
  expect_lte(max(trial_frame$recruitment_date), a)

  expect_lte(max(trial_frame$recruitment_date + trial_frame$Tstop), a + f)

  test_times <- pmin(trial_frame$Tstop, a + f - trial_frame$recruitment_date)
  test_times <- test_times[which(trial_frame$Tstart <= a + f - trial_frame$recruitment_date)]
  expect_equal(trial_frame$Tstop, test_times)

})
