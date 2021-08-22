test_that(".convolve_delay_distribution_vectors returns a vector whose elements sum up to 1", {
  N <- 100
  shapes <- stats::runif(2 * N, min = 0, max = 10)
  scales <- stats::runif(2 * N, min = 0, max = 10)

  distribution_list <- lapply(1:length(shapes), function(i) {
    return(list(name = "gamma", shape = shapes[i], scale = scales[i]))
  })

  delay_distribution_vectors <- lapply(distribution_list, function(x) {
    build_delay_distribution(x,
      max_quantile = 0.9999
    )
  })

  convolved_distribution_vectors <- sapply(1:N, function(x) {
    .convolve_delay_distribution_vectors(
      delay_distribution_vectors[[2 * x]],
      delay_distribution_vectors[[2 * x - 1]]
    )
  })

  max_difference_to_1 <- max(abs(sapply(delay_distribution_vectors, sum) - 1))

  expect_equal(max_difference_to_1, 0, tolerance = 1E-4)
})

test_that(".convolve_delay_distribution_vectors returns correct output on a simple example", {
  delay_distribution_vector_1 <- c(0, 0.25, 0.1, 0.65)
  delay_distribution_vector_2 <- c(0.2, 0.2, 0.3, 0.3)
  ref_convolved_output <- c(0, 0.05, 0.07, 0.225, 0.235, 0.225, 0.195, 0)

  convolved_output <- .convolve_delay_distribution_vectors(
    delay_distribution_vector_1,
    delay_distribution_vector_2
  )

  expect_equal(convolved_output, ref_convolved_output, tolerance = 1E-4)
})

test_that("convolve_delays returns same output as empirical method of convoluting gammas", {
  shape_incubation <- 3.2
  scale_incubation <- 1.3

  incubation_delay <- list(
    name = "gamma",
    shape = shape_incubation,
    scale = scale_incubation
  )

  shape_onset_to_report <- 2.7
  scale_onset_to_report <- 1.6

  onset_to_report_delay <- list(
    name = "gamma",
    shape = shape_onset_to_report,
    scale = scale_onset_to_report
  )

  convolved_output <- convolve_delays(
    delays = list(incubation_delay, onset_to_report_delay)
  )

  empirical_convolution_result <- c(
    0, 9e-04, 0.00947, 0.03214, 0.06438, 0.09523, 0.11566,
    0.12255, 0.11742, 0.1043, 0.08721, 0.06951, 0.05329,
    0.03951, 0.02841, 0.01991, 0.01371, 0.00925, 0.00616,
    0.00402, 0.0026, 0.00165, 0.00105, 0.00065, 0.00039,
    0.00025, 0.00015, 9e-05, 6e-05, 3e-05, 2e-05, 1e-05,
    1e-05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  )

  padded_convolved_output <- c(convolved_output, rep(0, times = length(empirical_convolution_result) - length(convolved_output)))

  absolute_diff <- abs(padded_convolved_output - empirical_convolution_result)

  expect_equal(max(absolute_diff), 0, tolerance = 1E-3)
})

test_that(".convolve_delay_distribution_vector_with_matrix returns correct output on a simple example", {
  vector_a <- c(0.2, 0.3, 0.5)
  matrix_b <- matrix(c(
    0.1, 0, 0,
    0.3, 0.2, 0,
    0.6, 0.4, 0.15
  ),
  nrow = 3,
  ncol = 3,
  byrow = TRUE
  )

  ref_convolved_matrix_vector_first <- matrix(c(
    0.02, 0, 0, 0, 0,
    0.09, 0.02, 0, 0, 0,
    0.26, 0.09, 0.02, 0, 0,
    0.33, 0.31, 0.12, 0.04, 0,
    0.30, 0.38, 0.315, 0.125, 0.03
  ),
  nrow = 5,
  ncol = 5,
  byrow = TRUE
  )

  ref_convolved_matrix_vector_last <- matrix(c(
    0.02, 0, 0,
    0.09, 0.04, 0,
    0.26, 0.14, 0.03
  ),
  nrow = 3,
  ncol = 3,
  byrow = TRUE
  )

  convolved_matrix_vector_first <- .convolve_delay_distribution_vector_with_matrix(
    vector_a = vector_a,
    matrix_b = matrix_b,
    vector_first = T
  )

  convolved_matrix_vector_last <- .convolve_delay_distribution_vector_with_matrix(
    vector_a = vector_a,
    matrix_b = matrix_b,
    vector_first = F
  )


  expect_equal(convolved_matrix_vector_first, ref_convolved_matrix_vector_first)
  expect_equal(convolved_matrix_vector_last, ref_convolved_matrix_vector_last)
})

test_that(".convolve_delay_distribution_vector_with_matrix returns valid output", {
  vector_a <- c(0.2, 0.3, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  vector_b <- c(0.3, 0.13, 0.42, 0.14, 0.01)
  matrix_b <- .get_matrix_from_single_delay_distr(vector_b, N = 20)

  convolved_matrix_vector_first <- .convolve_delay_distribution_vector_with_matrix(
    vector_a = vector_a,
    matrix_b = matrix_b,
    vector_first = T
  )

  convolved_matrix_vector_last <- .convolve_delay_distribution_vector_with_matrix(
    vector_a = vector_a,
    matrix_b = matrix_b,
    vector_first = F
  )

  expect_delay_matrix_sums_lte_1(convolved_matrix_vector_first, full_cols = 10)
  expect_delay_matrix_sums_lte_1(convolved_matrix_vector_last, full_cols = 10)
})

test_that(".convolve_delay_distribution_matrices returns valid output", {
  skip("Function is not ready yet.")

  vector_a <- c(0.21, 0.14, 0.17, 0.09, 0.01, 0.27, 0.11)
  vector_b <- c(0.3, 0.13, 0.42, 0.14, 0.01)
  matrix_a <- .get_matrix_from_single_delay_distr(vector_a, N = 30)
  matrix_b <- .get_matrix_from_single_delay_distr(vector_b, N = 30)

  convolved_matrix_ab <- .convolve_delay_distribution_matrices(
    matrix_a = matrix_a,
    matrix_b = matrix_b
  )

  convolved_matrix_ba <- .convolve_delay_distribution_matrices(
    matrix_a = matrix_b,
    matrix_b = matrix_a
  )

  expect_delay_matrix_sums_lte_1(convolved_matrix_ab, full_cols = 10)
  expect_delay_matrix_sums_lte_1(convolved_matrix_ba, full_cols = 10)
})

test_that(".convolve_delays is consistent on convolving several vectors and matrices", {
  vector_a <- c(0.2, 0.3, 0.5)
  matrix_b <- matrix(c(
    0.1, 0, 0,
    0.3, 0.2, 0,
    0.6, 0.4, 0.15
  ),
  nrow = 3,
  ncol = 3,
  byrow = TRUE
  )
  vector_c <-  c(0.2, 0.2, 0.3, 0.3)

  convolved_output <- convolve_delays(
    delays = list(vector_a, matrix_b, vector_c)
  )

  ref_result <- matrix(
    c(0.004, 0.000, 0.000, 0.000, 0.000,
  0.022, 0.004, 0.000, 0.000, 0.000,
  0.076, 0.022, 0.004, 0.000, 0.000,
  0.151, 0.086, 0.028, 0.008, 0.000,
  0.231, 0.171, 0.093, 0.033, 0.006),
  nrow = 5,
  ncol = 5,
  byrow = TRUE)

  expect_equal(convolved_output, ref_result, tolerance = 1E-3)
})


test_that(".convolve_delays can work with a single delay as input", {
  vector_a <- c(0.2, 0.3, 0.5)

  convolved_output <- convolve_delays(
    delays = vector_a
  )

  ref_result <- c(0.2, 0.3, 0.5)

  expect_equal(convolved_output, ref_result, tolerance = 1E-3)

  shape_incubation <- 3.2
  scale_incubation <- 1.3

  incubation_delay <- list(
    name = "gamma",
    shape = shape_incubation,
    scale = scale_incubation
  )

  convolved_output <- convolve_delays(
    delays = incubation_delay
  )

  ref_result <- c(0,0.08,0.17,0.2,0.17,0.13,0.09,0.06,0.04,0.02,0.01,0.01,0,0,0,0,0)

  expect_equal(convolved_output, ref_result, tolerance = 5E-2)
})
