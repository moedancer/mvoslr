library(mstate)
library(mvoslr)

test_that("execution for multiple sample sizes has the same outputs as simple evaluation", {

  my_seed <- 2

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

  msm_data_example_full <- data.frame(id = c(1,1,2,2,2,3,3,4,4,4), Tstart = c(0,0,0,0,0.5,0,0,0,0,1),
                                      Tstop = c(0.8,0.8,0.5,0.5,1.2,1.3,1.3,1,1,1.1),
                                      duration = c(0.8,0.8,0.5,0.5,0.7,1.3,1.3,1,1,0.1),
                                      from = c(1,1,1,1,2,1,1,1,1,2), to = c(2,3,2,3,3,2,3,2,3,3),
                                      status = c(0,1,1,0,1,0,0,1,0,0), trans = c(1,2,1,2,3,1,2,1,2,3),
                                      recruitment_date = c(0,0,0.3,0.3,0.3,0.7,0.7,0.9,0.9,0.9))

  msm_data_example_reduced <- msm_data_example_full[which(msm_data_example_full$id %in% 1:3), ]

  simple_result_full <- execution_mvoslr(msm_data_example_full, analysis_dates = analysis_dates_example, current_analysis = 2,
                                         reference_model = reference_model_example, events = events_example)
  simple_result_reduced <- execution_mvoslr(msm_data_example_reduced, analysis_dates = analysis_dates_example, current_analysis = 2,
                                            reference_model = reference_model_example, events = events_example)

  # For this test, we need to run 'execution_mvoslr_by_n' separately as subsets are chosen randomly
  enhanced_result_full <- execution_mvoslr_by_n(msm_data = msm_data_example_full, analysis_dates = analysis_dates_example, accrual_duration = 1,
                                                current_analysis = 2, reference_model = reference_model_example, events = events_example, sample_sizes = c(4))
  enhanced_result_reduced <- execution_mvoslr_by_n(msm_data = msm_data_example_reduced, analysis_dates = analysis_dates_example, accrual_duration = 1,
                                                current_analysis = 2, reference_model = reference_model_example, events = events_example, sample_sizes = c(3))

  expect_equal(unname(enhanced_result_full$raw_martingale[,,1]),
               unname(simple_result_full$raw_martingale))
  expect_equal(unname(enhanced_result_reduced$raw_martingale[,,1]),
               unname(simple_result_reduced$raw_martingale))

  expect_equal(unname(enhanced_result_reduced$stagewise_p_values[,1]),
               unname(simple_result_reduced$stagewise_p_values))
  expect_equal(unname(enhanced_result_full$stagewise_p_values[,1]),
               unname(simple_result_full$stagewise_p_values))

})

test_that("execution for multiple accrual durations has the same outputs as simple evaluation", {
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

  msm_data_example_full <- data.frame(id = c(1,1,2,2,2,3,3,4,4,4), Tstart = c(0,0,0,0,0.5,0,0,0,0,1),
                                      Tstop = c(0.8,0.8,0.5,0.5,1.2,1.3,1.3,1,1,1.1),
                                      duration = c(0.8,0.8,0.5,0.5,0.7,1.3,1.3,1,1,0.1),
                                      from = c(1,1,1,1,2,1,1,1,1,2), to = c(2,3,2,3,3,2,3,2,3,3),
                                      status = c(0,1,1,0,1,0,0,1,0,0), trans = c(1,2,1,2,3,1,2,1,2,3),
                                      recruitment_date = c(0,0,0.3,0.3,0.3,0.7,0.7,0.9,0.9,0.9))
  msm_data_example_reduced <- msm_data_example_full[which(msm_data_example_full$id <= 3), ]

  simple_result_full <- execution_mvoslr(msm_data_example_full, analysis_dates = analysis_dates_example, current_analysis = 2,
                                         reference_model = reference_model_example, events = events_example)
  simple_result_reduced <- execution_mvoslr(msm_data_example_reduced, analysis_dates = analysis_dates_example, current_analysis = 2,
                                            reference_model = reference_model_example, events = events_example)

  enhanced_result <- execution_mvoslr_by_a(msm_data = msm_data_example_full, analysis_dates = analysis_dates_example, accrual_durations = c(0.8, 1),
                                           current_analysis = 2, transition_matrix = tmat_example, cum_hazard_functions = cum_hazards_example,
                                           model_type = model_type_example, events = events_example)

  expect_equal(unname(enhanced_result$raw_martingale[,,1]),
               unname(simple_result_reduced$raw_martingale))
  expect_equal(unname(enhanced_result$raw_martingale[,,2]),
               unname(simple_result_full$raw_martingale))

  expect_equal(unname(enhanced_result$stagewise_p_values[,1]),
               unname(simple_result_reduced$stagewise_p_values))
  expect_equal(unname(enhanced_result$stagewise_p_values[,2]),
               unname(simple_result_full$stagewise_p_values))

})

test_that("execution for multiple accrual durations with fixed follow-up duration has the same outputs as simple evaluation", {
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
  accrual_durations_example <- c(0.8, 1)
  follow_up_example <- 1
  events_example <- list(c(2,3), c(3))
  names(events_example) <- c("PFS", "OS")

  msm_data_example_full <- data.frame(id = c(1,1,2,2,2,3,3,4,4,4), Tstart = c(0,0,0,0,0.5,0,0,0,0,1),
                                      Tstop = c(0.8,0.8,0.5,0.5,1.2,1.3,1.3,1,1,1.1),
                                      duration = c(0.8,0.8,0.5,0.5,0.7,1.3,1.3,1,1,0.1),
                                      from = c(1,1,1,1,2,1,1,1,1,2), to = c(2,3,2,3,3,2,3,2,3,3),
                                      status = c(0,1,1,0,1,0,0,1,0,0), trans = c(1,2,1,2,3,1,2,1,2,3),
                                      recruitment_date = c(0,0,0.3,0.3,0.3,0.7,0.7,0.9,0.9,0.9))
  msm_data_example_reduced <- data.frame(id = c(1,1,2,2,2,3,3), Tstart = c(0,0,0,0,0.5,0,0),
                                         Tstop = c(0.8,0.8,0.5,0.5,1.2,1.1,1.1),
                                         duration = c(0.8,0.8,0.5,0.5,0.7,1.1,1.1),
                                         from = c(1,1,1,1,2,1,1), to = c(2,3,2,3,3,2,3),
                                         status = c(0,1,1,0,1,0,0), trans = c(1,2,1,2,3,1,2),
                                         recruitment_date = c(0,0,0.3,0.3,0.3,0.7,0.7))

  simple_result_full <- execution_mvoslr(msm_data_example_full, analysis_dates = analysis_dates_example, current_analysis = 2,
                                         reference_model = reference_model_example, events = events_example)
  simple_result_reduced <- execution_mvoslr(msm_data_example_reduced, analysis_dates = analysis_dates_example, current_analysis = 2,
                                            reference_model = reference_model_example, events = events_example)

  enhanced_result <- execution_mvoslr_fixed_fu(msm_data = msm_data_example_full, interim_analysis_dates = head(analysis_dates_example, -1),
                                               accrual_durations = accrual_durations_example,
                                               follow_up = follow_up_example, current_analysis = 2, transition_matrix = tmat_example,
                                               cum_hazard_functions = cum_hazards_example,
                                               model_type = model_type_example, events = events_example)

  expect_equal(unname(enhanced_result$raw_martingale[,,1]),
               unname(simple_result_reduced$raw_martingale))
  expect_equal(unname(enhanced_result$raw_martingale[,,2]),
               unname(simple_result_full$raw_martingale))

  expect_equal(unname(enhanced_result$stagewise_p_values[,1]),
               unname(simple_result_reduced$stagewise_p_values))
  expect_equal(unname(enhanced_result$stagewise_p_values[,2]),
               unname(simple_result_full$stagewise_p_values))

})
