#TODO test LOESS function
# 1) always positive values
# 2) normalized total
# 3) constant output if constant input
# 4) add a "reference" test for a more complicated use-case
# 5) test that min and max of smoothed values are bounded by min and max of noisy values

#TODO continue

smoothing_methods_tested <- c("LOESS")

test_that("smooth_incidence output values are positive ", {
  random_values <- sample.int(1000, size = 100, replace= TRUE)

  sapply(smoothing_methods_tested, function(x) {
        smoothed_values <- .get_values(smooth_incidence(random_values, smoothing_method = x))
         expect_gte(min(smoothed_values), 0)})
})

test_that("smooth_incidence output values are constant if input is constant ", {
  constant_values <- rep(5.3, times = 100)

  sapply(smoothing_methods_tested, function(x) {
    smoothed_values <- .get_values(smooth_incidence(constant_values, smoothing_method = x))
    expect_equal(smoothed_values, constant_values)})
})

test_that("smooth_incidence output values sum up to the total of the original incidence ", {
  random_values <- sample.int(1000, size = 100, replace= TRUE)

  sapply(smoothing_methods_tested, function(x) {
    smoothed_values <- .get_values(smooth_incidence(random_values, smoothing_method = x))
    expect_equal(sum(smoothed_values), sum(random_values))})
})
