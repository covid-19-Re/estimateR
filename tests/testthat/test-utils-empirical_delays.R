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
  skip_on_cran()
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
  expect_equal(max(result$scale_rmse, 0.3808837), 0.3808837, tolerance = 1E-2) # kept the lower time_steps value to reduce running time


  # Testing delay matrix with data sampled from a different gamma distribution for each timestep
  original_distribution_shapes <- sample(seq(3.9, 7.1, by = 0.1), time_steps, replace = TRUE)
  original_distribution_scales <- sample(seq(2.9, 6.1, by = 0.1), time_steps, replace = TRUE)
  set.seed(1)
  result <- .delay_distribution_matrix_rmse_compute(original_distribution_shapes, original_distribution_scales, nr_distribution_samples)
  expect_equal(max(result$shape_rmse, 0.2932824), 0.2932824, tolerance = 1E-2)
  expect_equal(max(result$scale_rmse, 0.2883487), 0.2883487, tolerance = 1E-2)
})
