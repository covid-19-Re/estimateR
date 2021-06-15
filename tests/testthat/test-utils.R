#TODO test that merge outputs works with inputs of

test_that(".get_module_input() deals with well-formatted input", {
  toy_incidence <- c(1,2,1,0)
  data_1 <- list(values = toy_incidence, index_offset = -4)

  expect_identical(.get_module_input(toy_incidence), list(values = toy_incidence, index_offset = 0))
  expect_identical(.get_module_input(data_1), data_1)
})

test_that(".get_module_input() checks input format", {
  toy_incidence <- c(1,23,1,50)
  data_1 <- list(values = toy_incidence, index_offset = -2, baz = 6)
  invalid_input_1 <- list(a = toy_incidence, b = -2, c = 6)

  expect_identical(.get_module_input(data_1), list(values = toy_incidence, index_offset = -2))
  expect_error(.get_module_input(invalid_input_1))
})

test_that(".get_module_output() deals with well-formatted input", {
  results <- c(1,2,3,-2,0,0)
  input_1 <- list(values = c(1,1,1,1,1), index_offset = 0)
  input_2 <- list(values = c(1,1,1,1,1), index_offset = -1)

  expect_identical(.get_module_output(results, input_1), list(values = results, index_offset = 0))
  expect_identical(.get_module_output(results, input_2), list(values = results, index_offset = -1))
  expect_identical(.get_module_output(results, input_2, offset = 3), list(values = results, index_offset = 2))
})

test_that(".get_module_output() checks input format", {
  #TODO make test
})

test_that("correct_for_partially_observed_data() is correct on a simple example", {
  toy_incidence <- c(9,3,4,7,8,0,2,1,0)
  toy_delay <- c(0.05,0.05,0.3,0.2,0.1,0.3)

  corrected_result <- correct_for_partially_observed_data(toy_incidence, toy_delay, cutoff_observation_probability = 0.11)

  ref_values <- c(9,3,4,7,11.429,0,5)
  expect_equal(corrected_result$values, ref_values, tolerance = 0.01)
})
