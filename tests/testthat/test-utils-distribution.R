# TODO TEST that:
# 1) build_delay_distribution throws error when unsupported distribution_type is thrown in
# and when unsuitable parameter values are thrown in (not numeric, or negative values for instance)

test_that("build_delay_distribution returns a vector whose elements sum up to 1", {
  N <- 100
  shapes <- stats::runif(N, min = 0, max = 10)
  scales <- stats::runif(N, min = 0, max = 10)

  distribution_list <- lapply(1:length(shapes), function(i) {
    return(list(name = "gamma", shape = shapes[i], scale = scales[i]))
  })

  delay_distribution_vectors <- lapply(distribution_list, function(x) {
    build_delay_distribution(x,
      max_quantile = 0.9999
    )
  })

  max_difference_to_1 <- max(abs(sapply(delay_distribution_vectors, sum) - 1))

  expect_equal(max_difference_to_1, 0, tolerance = 1E-4)
})


test_that(".get_delay_matrix_from_delay_distributions returns valid output", {
  N <- 100

  shapes <- stats::runif(N, min = 0, max = 10)
  scales <- stats::runif(N, min = 0, max = 10)

  distribution_list <- lapply(1:length(shapes), function(i) {
    return(list(name = "gamma", shape = shapes[i], scale = scales[i]))
  })

  matrix_result <- .get_delay_matrix_from_delay_distributions(
    distributions = distribution_list,
    max_quantile = 0.999
  )

  # Check that all columns sum up to less than one.
  expect_delay_matrix_sums_lte_1(matrix_result, full_cols = 0)
})
