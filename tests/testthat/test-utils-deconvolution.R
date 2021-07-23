# TODO reorganize: check that all utilities tested are still in utils-deconvolution.R (create additional testing files otherwise and move tests there)

# TODO TEST that:
# 1) build_delay_distribution throws error when unsupported distribution_type is thrown in
# and when unsuitable parameter values are thrown in (not numeric, or negative values for instance)
# 2) get_matrix_empirical_waiting_time_distr: last few columns are constant when they should be (on a simple example)

expect_delay_matrix_sums_lte_1 <- function(matrix, full_cols = 0, tolerance = 1E-3) {
  if (full_cols > 0 && full_cols < ncol(matrix)) {
    sums_full_cols <- apply(matrix[, 1:full_cols], MARGIN = 2, FUN = sum)
    expect_equal(sums_full_cols, rep(1, times = length(sums_full_cols)), tolerance = tolerance)
  }

  sum_all_cols <- apply(matrix, MARGIN = 2, FUN = sum)
  expect_lte(max(abs(sum_all_cols)), 1)
}

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

test_that("convolve_delay_inputs returns same output as empirical method of convoluting gammas", {
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

  # TODO fix this (need to pass n_report_time_steps)
  convolved_output <- convolve_delay_inputs(
    delay_incubation = incubation_delay,
    delay_onset_to_report = onset_to_report_delay
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

test_that(".get_delay_matrix_from_delay_distribution_parms returns valid output", {
  N <- 100

  shapes <- stats::runif(N, min = 0, max = 10)
  scales <- stats::runif(N, min = 0, max = 10)

  distribution_list <- lapply(1:length(shapes), function(i) {
    return(list(name = "gamma", shape = shapes[i], scale = scales[i]))
  })

  matrix_result <- .get_delay_matrix_from_delay_distribution_parms(
    distributions = distribution_list,
    max_quantile = 0.999
  )

  # Check that all columns sum up to less than one.
  expect_delay_matrix_sums_lte_1(matrix_result, full_cols = 0)
})

test_that(".get_matrix_from_empirical_delay_distr returns valid output", {
  # First toy data test
  ref_date <- as.Date("2020-03-01")
  time_series_length <- 100

  report_delays <- sample(c(0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 7, 8, 9, 10), 1000, replace = T)
  event_dates <- sample(seq.Date(from = ref_date, length.out = time_series_length, by = "day"), 1000, replace = T)

  empirical_delay_data <- tibble::tibble(
    event_date = event_dates,
    report_delay = report_delays
  ) %>%
    dplyr::arrange(event_date)

  empirical_matrix <- get_matrix_from_empirical_delay_distr(
    empirical_delays = empirical_delay_data,
    ref_date = ref_date,
    n_report_time_steps = 90,
    time_step = "day",
    min_number_cases = 10,
    upper_quantile_threshold = 0.99,
    fit = "none"
  )

  expect_delay_matrix_sums_lte_1(empirical_matrix, full_cols = 50)

  # Second toy data test
  ref_date <- as.Date("2020-04-01")
  n_days <- 50
  delay_increase <- 1.5
  shape_initial_delay <- 6
  scale_initial_delay <- 1.5
  distribution_initial_delay <- list(name = "gamma", shape = shape_initial_delay, scale = scale_initial_delay)
  seed <- 734


  generated_empirical_delays <- .generate_delay_data(
    origin_date = ref_date,
    n_time_steps = n_days,
    ratio_delay_end_to_start = 1.5,
    distribution_initial_delay = distribution_initial_delay,
    seed = seed
  )

  empirical_delays_matrix <- get_matrix_from_empirical_delay_distr(
    empirical_delays = generated_empirical_delays,
    ref_date = ref_date,
    n_report_time_steps = 50,
    fit = "none",
    min_number_cases = 5
  )


  expect_delay_matrix_sums_lte_1(empirical_delays_matrix, full_cols = 20)
})

test_that(".get_matrix_from_empirical_delay_distr handles returning data over a full number of weeks", {
  # First toy data test
  ref_date <- as.Date("2020-03-01")
  time_series_length <- 100

  report_delays <- sample(c(0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 7, 8, 9, 10), 1000, replace = T)
  event_dates <- sample(seq.Date(from = ref_date, length.out = time_series_length, by = "day"), 1000, replace = T)

  empirical_delay_data <- tibble::tibble(
    event_date = event_dates,
    report_delay = report_delays
  ) %>%
    dplyr::arrange(event_date)

  empirical_matrix <- get_matrix_from_empirical_delay_distr(
    empirical_delays = empirical_delay_data,
    ref_date = ref_date,
    n_report_time_steps = 90,
    time_step = "day",
    min_number_cases = 10,
    upper_quantile_threshold = 0.99,
    fit = "none",
    num_steps_in_a_unit = 7
  )

  expect_delay_matrix_sums_lte_1(empirical_matrix, full_cols = 50)
})

test_that(".get_matrix_from_empirical_delay_distr returns a matrix with the expected distributions when using fit = gamma", {
  nr_distribution_samples <- 500
  time_steps <- 30

  # Testing delay matrix with data sampled from constant gamma distribution
  original_distribution_shapes <- rep(6, time_steps)
  original_distribution_scales <- rep(5, time_steps)
  set.seed(1)
  result <- .delay_distribution_matrix_rmse_compute(original_distribution_shapes, original_distribution_scales, nr_distribution_samples)
  expect_equal(max(result$shape_rmse, 0.07829915), 0.07829915, tolerance = 1E-2)
  expect_equal(max(result$scale_rmse, 0.05670633), 0.05670633, tolerance = 1E-2)


  # Testing delay matrix with data sampled from two different gamma distributions
  original_distribution_shapes <- c(rep(3.5, time_steps / 2), rep(6.5, time_steps / 2))
  original_distribution_scales <- c(rep(2, time_steps / 2), rep(3, time_steps / 2))
  set.seed(1)
  result <- .delay_distribution_matrix_rmse_compute(original_distribution_shapes, original_distribution_scales, nr_distribution_samples)
  expect_equal(max(result$shape_rmse, 0.3193425), 0.3193425, tolerance = 1E-2) # the RMSE gets lower with more time_steps;
  expect_equal(max(result$scale_rmse, 0.3491427), 0.3491427, tolerance = 1E-2) # kept the lower time_steps value to reduce running time


  # Testing delay matrix with data sampled from a different gamma distribution for each timestep
  original_distribution_shapes <- sample(seq(3.9, 7.1, by = 0.1), time_steps, replace = TRUE)
  original_distribution_scales <- sample(seq(2.9, 6.1, by = 0.1), time_steps, replace = TRUE)
  set.seed(1)
  result <- .delay_distribution_matrix_rmse_compute(original_distribution_shapes, original_distribution_scales, nr_distribution_samples)
  expect_equal(max(result$shape_rmse, 0.1604508), 0.1604508, tolerance = 1E-2)
  expect_equal(max(result$scale_rmse, 0.2854936), 0.2854936, tolerance = 1E-2)
})
