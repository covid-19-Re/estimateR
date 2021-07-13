# TODO test LOESS function
# 1) always positive values
# 2) normalized total
# 3) constant output if constant input
# 4) add a "reference" test for a more complicated use-case
# 5) test that min and max of smoothed values are bounded by min and max of noisy values

smoothing_methods_tested <- c("LOESS")

test_that("smooth_incidence output values are positive ", {
  random_values <- sample.int(1000, size = 100, replace = TRUE)

  sapply(smoothing_methods_tested, function(x) {
    smoothed_values <- .get_values(smooth_incidence(random_values, smoothing_method = x))
    expect_gte(min(smoothed_values), 0)
  })
})

test_that("smooth_incidence output values sum up to the total of the original incidence ", {
  random_values <- sample.int(1000, size = 100, replace = TRUE)

  sapply(smoothing_methods_tested, function(x) {
    smoothed_values <- .get_values(smooth_incidence(random_values, smoothing_method = x))
    expect_equal(sum(smoothed_values), sum(random_values))
  })
})

test_that("smooth_incidence output values are constant if input is constant ", {
  # TODO unskip test when added way to deal with left-truncated incidence input
  skip("Left-truncated incidence input cannot be dealt with yet.")
  constant_values <- rep(5.3, times = 100)

  sapply(smoothing_methods_tested, function(x) {
    smoothed_values <- .get_values(smooth_incidence(constant_values, smoothing_method = x))
    expect_equal(smoothed_values, constant_values)
  })
})

test_that("smooth_incidence output values are bounded by bounds of noisy values", {
  random_values <- sample.int(1000, size = 100, replace = TRUE)

  sapply(smoothing_methods_tested, function(x) {
    smoothed_values <- .get_values(smooth_incidence(random_values, smoothing_method = x))
    expect_gte(min(smoothed_values), min(random_values))
    expect_lte(max(smoothed_values), max(random_values))
  })
})

test_that("smooth_incidence output stays consistent for LOESS method", {
  noisy_values <- c(272, 78, 859, 642, 411, 612, 192, 262, 399, 371, 69, 80, 221, 945, 198, 896, 705, 155, 498, 795)
  ref_smoothed_values <- c(213.481, 238.774, 265.581, 293.119, 317.855, 341.012, 365.157, 388.764, 410.31, 430.617, 451.25, 471.793, 491.829, 511.203, 530.235, 549.217, 568.446, 587.809, 607.097, 626.452)
  smoothed_values <- .get_values(smooth_incidence(noisy_values, smoothing_method = "LOESS"))
  expect_equal(smoothed_values, ref_smoothed_values, tolerance = 1E-3)
})
