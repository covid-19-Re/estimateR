test_that("correct_for_partially_observed_data() is correct on a simple example", {
  toy_incidence <- c(9, 3, 4, 7, 8, 0, 2, 1, 0)
  toy_delay <- c(0.05, 0.05, 0.3, 0.2, 0.1, 0.3)

  corrected_result <- correct_for_partially_observed_data(toy_incidence, toy_delay, cutoff_observation_probability = 0.11)

  ref_values <- c(9, 3, 4, 7, 11.429, 0, 5)
  expect_equal(corrected_result$values, ref_values, tolerance = 0.01)
})
