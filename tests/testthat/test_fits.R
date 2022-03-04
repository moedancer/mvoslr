library(mvoslr)
library(mstate)

test_that("parameters in time-homogeneous Markov model are correctly estimated", {
  test_data <- data.frame(id = 1:10, Tstart = rep(0,10), Tstop = 1:10,
                          duration = 1:10, from = rep(1,10), to = rep(2,10),
                          status = c(rep(0,5), rep(1,5)), trans = rep(1,10))
  test_trans <- trans.comprisk(1)
  expect_equal(unname(fit_thmm(test_data, test_trans)$parameters[2, 1]),
               unname(sum(test_data$duration)/sum(test_data$status)))
})


test_that("estimates for time-inhomogeneous and Semi-Markov model
          are the same in case of two-state model", {
  test_data <- data.frame(id = 1:10, Tstart = rep(0,10), Tstop = 1:10,
                          duration = 1:10, from = rep(1,10), to = rep(2,10),
                          status = c(rep(0,5), rep(1,5)), trans = rep(1,10))
  test_trans <- trans.comprisk(1)
  expect_equal(unname(fit_tihmm(test_data, test_trans)$parameters),
               unname(fit_smm(test_data, test_trans)$parameters))
})

test_that("estimates for time-inhomogeneous and Semi-Markov model
          are the same in case of competing risks", {
            test_data <- data.frame(id = rep(1:9, each = 3),
                                    Tstart = rep(0, 27), Tstop = rep(1:9, each = 3),
                                    duration = rep(1:9, each = 3), from = rep(1, 27),
                                    to = rep(2:4, 9),
                                    status = rep(c(1,0,0,0,1,0,0,0,1), 3),
                                    trans = rep(1:3, 9))
            test_trans <- trans.comprisk(3)
            expect_equal(unname(fit_tihmm(test_data, test_trans)$parameters),
                         unname(fit_smm(test_data, test_trans)$parameters))
})
