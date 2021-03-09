#TODO test LOESS function
# 1) always positive values
# 2) normalized total
# 3) constant output if constant input
# 4) add a "reference" test for a more complicated use-case
# 5) test that min and max of smoothed values are bounded by min and max of noisy values

#TODO continue

smoothing_methods_tested <- c("LOESS")

test_that("smooth_incidence output values is positive  ", {
  expect_equal(2 * 2, 4)
})
