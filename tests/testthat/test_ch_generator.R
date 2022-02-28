library(mvoslr)

test_that("cumulative hazard functions with exponential shape are generated properly", {
  my_rate <- 10
  my_t <- 5
  expect_equal(get_cum_haz_fct_exp(my_rate)(my_t),
               my_rate * my_t)
})

test_that("cumulative hazard functions with Weibull shape are generated properly", {
  my_shape <- 2
  my_scale <- 10
  my_t <- 3
  expect_equal(get_cum_haz_fct_weibull(c(my_shape, my_scale))(my_t),
               -log(1 - pweibull(my_t, shape = my_shape, scale = my_scale)))
})
