library(mvoslr)
library(mstate)

test_that("power function correctly aggregates results from single trials", {

  my_seed <- 101

  tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
  number_of_trans_example <- dim(to.trans2(tmat_example))[1]
  model_type_example <- "SM"
  cumhaz_12_example <- function(t) t^1.1
  cumhaz_13_example <- function(t) t^1.2
  cumhaz_23_example <- function(t) t^0.9
  cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
  analysis_dates_example <- c(1, 2)
  events_example <- list(c(2,3), c(3))
  names(events_example) <- c("PFS", "OS")
  accrual_duration_example <- 1
  sample_size_example <- 100
  #In this example, the alternative is specified via separate hazard ratios for each transition
  hazard_ratios_example <- c(1.4, 1.2, 1.35)

  cum_hazards_alternative_example <- vector(mode = "list", length = number_of_trans_example)
  for(trans in 1:number_of_trans_example){
    cum_hazards_alternative_example[[trans]] <- function(x){
      (1/hazard_ratios_example[trans]) * do.call(cum_hazards_example[[trans]], list(x))
    }
  }
  cumhaz_alternative <- discretize_functions(cum_hazards_alternative_example,
                                             analysis_dates_example[2], time_steps = 100)

  set.seed(my_seed)

  sim_frame_1 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_1 <- msm_to_trial_data(sim_frame_1, accrual_duration_example,
                                   analysis_dates_example[2] - accrual_duration_example)
  result_1 <- execution_mvoslr(msm_data = sim_frame_1, analysis_dates = analysis_dates_example,
                               current_analysis = 2, transition_matrix = tmat_example,
                               cum_hazard_functions = cum_hazards_example,
                               model_type = model_type_example, events = events_example)
  decision_1 <- !is.na(result_1$rejection_stage)

  sim_frame_2 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_2 <- msm_to_trial_data(sim_frame_2, accrual_duration_example,
                                   analysis_dates_example[2] - accrual_duration_example)
  result_2 <- execution_mvoslr(msm_data = sim_frame_2, analysis_dates = analysis_dates_example,
                               current_analysis = 2, transition_matrix = tmat_example,
                               cum_hazard_functions = cum_hazards_example,
                               model_type = model_type_example, events = events_example)
  decision_2 <- !is.na(result_2$rejection_stage)

  sim_frame_3 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_3 <- msm_to_trial_data(sim_frame_3, accrual_duration_example,
                                   analysis_dates_example[2] - accrual_duration_example)
  result_3 <- execution_mvoslr(msm_data = sim_frame_3, analysis_dates = analysis_dates_example,
                               current_analysis = 2, transition_matrix = tmat_example,
                               cum_hazard_functions = cum_hazards_example,
                               model_type = model_type_example, events = events_example)
  decision_3 <- !is.na(result_3$rejection_stage)

  set.seed(my_seed)

  power_result <- power_mvoslr(transition_matrix = tmat_example, model_type = model_type_example,
                               events = events_example, cum_hazard_functions_h0 = cum_hazards_example,
                               analysis_dates = analysis_dates_example,
                               accrual_duration = accrual_duration_example, sample_size = sample_size_example,
                               hazard_ratios = hazard_ratios_example, simulation_runs = 3)

  expect_equal(unname((result_1$raw_martingale + result_2$raw_martingale + result_3$raw_martingale)/3),
               power_result$means)
  expect_equal(mean(c(decision_1, decision_2, decision_3)),
               power_result$power)

})

test_that("power function with variable accrual duration correctly aggregates results from single trials", {

  my_seed <- 102

  tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
  number_of_trans_example <- dim(to.trans2(tmat_example))[1]
  model_type_example <- "SM"
  cumhaz_12_example <- function(t) t^1.1
  cumhaz_13_example <- function(t) t^1.2
  cumhaz_23_example <- function(t) t^0.9
  cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
  analysis_dates_example <- c(1, 2)
  events_example <- list(c(2,3), c(3))
  names(events_example) <- c("PFS", "OS")
  accrual_durations_example <- c(0.5, 1)
  max_accrual_duration_example <- max(accrual_durations_example)
  sample_size_example <- 100
  recruitment_speed_example <- 100
  #In this example, the alternative is specified via separate hazard ratios for each transition
  hazard_ratios_example <- c(1.4, 1.2, 1.35)

  cum_hazards_alternative_example <- vector(mode = "list", length = number_of_trans_example)
  for(trans in 1:number_of_trans_example){
    cum_hazards_alternative_example[[trans]] <- function(x){
      (1/hazard_ratios_example[trans]) * do.call(cum_hazards_example[[trans]], list(x))
    }
  }
  cumhaz_alternative <- discretize_functions(cum_hazards_alternative_example,
                                             analysis_dates_example[2], time_steps = 100)

  set.seed(my_seed)

  sim_frame_1 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_1 <- msm_to_trial_data(sim_frame_1, accrual_durations_example[2],
                                   analysis_dates_example[2] - accrual_durations_example[2])

  result_1 <- execution_mvoslr_by_a(msm_data = sim_frame_1, analysis_dates = analysis_dates_example,
                                    current_analysis = 2, transition_matrix = tmat_example,
                                    cum_hazard_functions = cum_hazards_example,
                                    model_type = model_type_example, events = events_example,
                                    accrual_durations = accrual_durations_example)

  decision_1_long <- !is.na(result_1$rejection_stage[1])
  decision_1_short <- !is.na(result_1$rejection_stage[2])

  sim_frame_2 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_2 <- msm_to_trial_data(sim_frame_2, accrual_durations_example[2],
                                   analysis_dates_example[2] - accrual_durations_example[2])

  result_2 <- execution_mvoslr_by_a(msm_data = sim_frame_2, analysis_dates = analysis_dates_example,
                                    current_analysis = 2, transition_matrix = tmat_example,
                                    cum_hazard_functions = cum_hazards_example,
                                    model_type = model_type_example, events = events_example,
                                    accrual_durations = accrual_durations_example)

  decision_2_long <- !is.na(result_2$rejection_stage[1])
  decision_2_short <- !is.na(result_2$rejection_stage[2])

  sim_frame_3 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_3 <- msm_to_trial_data(sim_frame_3, accrual_durations_example[2],
                                   analysis_dates_example[2] - accrual_durations_example[2])

  result_3 <- execution_mvoslr_by_a(msm_data = sim_frame_3, analysis_dates = analysis_dates_example,
                                    current_analysis = 2, transition_matrix = tmat_example,
                                    cum_hazard_functions = cum_hazards_example,
                                    model_type = model_type_example, events = events_example,
                                    accrual_durations = accrual_durations_example)

  decision_3_long <- !is.na(result_3$rejection_stage[1])
  decision_3_short <- !is.na(result_3$rejection_stage[2])

  set.seed(my_seed)
  (result_1$raw_martingale[,,1] + result_2$raw_martingale[,,1] + result_3$raw_martingale[,,1])/3
  (result_1$raw_martingale[,,2] + result_2$raw_martingale[,,2] + result_3$raw_martingale[,,2])/3

  power_result <- power_mvoslr_by_a(transition_matrix = tmat_example, model_type = model_type_example,
                                    events = events_example, cum_hazard_functions_h0 = cum_hazards_example,
                                    analysis_dates = analysis_dates_example,
                                    accrual_durations = accrual_durations_example, recruitment_speed = recruitment_speed_example,
                                    hazard_ratios = hazard_ratios_example, simulation_runs = 3)

  expect_equal(unname((result_1$raw_martingale[,,1] + result_2$raw_martingale[,,1] + result_3$raw_martingale[,,1])/3),
               unname(power_result$means[,,1]))
  expect_equal(unname((result_1$raw_martingale[,,2] + result_2$raw_martingale[,,2] + result_3$raw_martingale[,,2])/3),
               unname(power_result$means[,,2]))
  expect_equal(unname(mean(c(decision_1_short, decision_2_short, decision_3_short))),
               unname(power_result$power[1]))
  expect_equal(unname(mean(c(decision_1_long, decision_2_long, decision_3_long))),
               unname(power_result$power[2]))

})

test_that("power function with variable accrual duration and fixed follow-up duration
          correctly aggregates results from single trials", {

  my_seed <- 103

  tmat_example <- transMat(x = list(c(2,3),c(3),c()), names = c("a", "b", "c"))
  number_of_trans_example <- dim(to.trans2(tmat_example))[1]
  model_type_example <- "M"
  cumhaz_12_example <- function(t) t^1
  cumhaz_13_example <- function(t) t^2
  cumhaz_23_example <- function(t) t^3
  cum_hazards_example <- list(cumhaz_12_example, cumhaz_13_example, cumhaz_23_example)
  interim_analysis_dates_example <- c(0.5, 1)
  events_example <- list(c(2,3), c(3))
  names(events_example) <- c("PFS", "OS")
  accrual_durations_example <- c(0.5, 1)
  follow_up_example <- 1
  max_accrual_duration_example <- max(accrual_durations_example)
  sample_size_example <- 100
  recruitment_speed_example <- 100
  #In this example, the alternative is specified via separate hazard ratios for each transition
  hazard_ratios_example <- c(1.1, 1.5, 0.8)

  cum_hazards_alternative_example <- vector(mode = "list", length = number_of_trans_example)
  for(trans in 1:number_of_trans_example){
    cum_hazards_alternative_example[[trans]] <- function(x){
      (1/hazard_ratios_example[trans]) * do.call(cum_hazards_example[[trans]], list(x))
    }
  }
  cumhaz_alternative <- discretize_functions(cum_hazards_alternative_example,
                                             max(accrual_durations_example) + follow_up_example, time_steps = 100)

  set.seed(my_seed)

  sim_frame_1 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_1 <- msm_to_trial_data(sim_frame_1, accrual_durations_example[2],
                                   follow_up_example + accrual_durations_example[2])

  result_1 <- execution_mvoslr_fixed_fu(msm_data = sim_frame_1, interim_analysis_dates = interim_analysis_dates_example,
                                        current_analysis = 3, transition_matrix = tmat_example, follow_up = follow_up_example,
                                        cum_hazard_functions = cum_hazards_example,
                                        model_type = model_type_example, events = events_example,
                                        accrual_durations = accrual_durations_example)

  decision_1_long <- !is.na(result_1$rejection_stage[1])
  decision_1_short <- !is.na(result_1$rejection_stage[2])

  sim_frame_2 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_2 <- msm_to_trial_data(sim_frame_2, accrual_durations_example[2],
                                   follow_up_example + accrual_durations_example[2])

  result_2 <- execution_mvoslr_fixed_fu(msm_data = sim_frame_2, interim_analysis_dates = interim_analysis_dates_example,
                                        current_analysis = 3, transition_matrix = tmat_example, follow_up = follow_up_example,
                                        cum_hazard_functions = cum_hazards_example,
                                        model_type = model_type_example, events = events_example,
                                        accrual_durations = accrual_durations_example)

  decision_2_long <- !is.na(result_2$rejection_stage[1])
  decision_2_short <- !is.na(result_2$rejection_stage[2])

  sim_frame_3 <- simulate_msm(tmat_example, model_type_example, cumhaz_alternative, sample_size_example)
  sim_frame_3 <- msm_to_trial_data(sim_frame_3, accrual_durations_example[2],
                                   follow_up_example + accrual_durations_example[2])

  result_3 <- execution_mvoslr_fixed_fu(msm_data = sim_frame_3, interim_analysis_dates = interim_analysis_dates_example,
                                        current_analysis = 3, transition_matrix = tmat_example, follow_up = follow_up_example,
                                        cum_hazard_functions = cum_hazards_example,
                                        model_type = model_type_example, events = events_example,
                                        accrual_durations = accrual_durations_example)

  decision_3_long <- !is.na(result_3$rejection_stage[1])
  decision_3_short <- !is.na(result_3$rejection_stage[2])

  set.seed(my_seed)
  (result_1$raw_martingale[,,1] + result_2$raw_martingale[,,1] + result_3$raw_martingale[,,1])/3
  (result_1$raw_martingale[,,2] + result_2$raw_martingale[,,2] + result_3$raw_martingale[,,2])/3

  power_result <- power_mvoslr_fixed_fu(transition_matrix = tmat_example, model_type = model_type_example,
                                        events = events_example, cum_hazard_functions_h0 = cum_hazards_example,
                                        interim_analysis_dates = interim_analysis_dates_example, follow_up = follow_up_example,
                                        accrual_durations = accrual_durations_example, recruitment_speed = recruitment_speed_example,
                                        hazard_ratios = hazard_ratios_example, simulation_runs = 3)

  expect_equal(unname((result_1$raw_martingale[,,1] + result_2$raw_martingale[,,1] + result_3$raw_martingale[,,1])/3),
               unname(power_result$means[,,1]))
  expect_equal(unname((result_1$raw_martingale[,,2] + result_2$raw_martingale[,,2] + result_3$raw_martingale[,,2])/3),
               unname(power_result$means[,,2]))
  expect_equal(unname(mean(c(decision_1_short, decision_2_short, decision_3_short))),
               unname(power_result$power[1]))
  expect_equal(unname(mean(c(decision_1_long, decision_2_long, decision_3_long))),
               unname(power_result$power[2]))

})
