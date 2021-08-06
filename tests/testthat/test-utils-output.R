test_that(".get_module_output() deals with well-formatted input", {
  results <- c(1, 2, 3, -2, 0, 0)
  input_1 <- list(values = c(1, 1, 1, 1, 1), index_offset = 0)
  input_2 <- list(values = c(1, 1, 1, 1, 1), index_offset = -1)

  expect_identical(.get_module_output(results, input_1), list(values = results, index_offset = 0))
  expect_identical(.get_module_output(results, input_2), list(values = results, index_offset = -1))
  expect_identical(.get_module_output(results, input_2, offset = 3), list(values = results, index_offset = 2))
})

test_that(".get_module_output() checks input format", {
  # TODO make test
})
