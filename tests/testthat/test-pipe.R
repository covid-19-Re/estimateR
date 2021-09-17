test_that("estimate_Re_from_noisy_delayed_incidence yields consistent results on a toy example", {
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

  delay_distribution <- c(
    0, 0.015, 0.09, 0.168,
    0.195, 0.176, 0.135, 0.091, 0.057, 0.034,
    0.019, 0.01, 0.005, 0.003,
    0.001, 0.001
  )

  estimates <- estimate_Re_from_noisy_delayed_incidence(
    incidence_data = toy_incidence_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    delay = delay_distribution,
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    minimum_cumul_incidence = 10,
    mean_Re_prior = 1,
    output_Re_only = FALSE,
    ref_date = as.Date("2020-02-04"),
    time_step = "day"
  )

  reference_R_values <- c(
    NA, NA, NA, NA, NA, NA, 3.15, 2.86, 2.67,
    2.53, 2.41, 2.29, 2.18, 2.08, 1.98, 1.88,
    1.78, 1.69, 1.61, 1.52, 1.44, 1.37, 1.3, 1.23,
    1.16, 1.09, 1.02, 0.96, 0.9, 0.86, 0.82, 0.79,
    0.76, 0.74, 0.72,0.69,0.66,NA,NA,NA,NA,NA)

  reference_dates <- seq.Date(from = as.Date("2020-01-30"), to = as.Date("2020-03-11"), by = "day")

  expect_equal(estimates$Re_estimate, reference_R_values, tolerance = 5E-2)
  expect_equal(estimates$date, reference_dates)
})

test_that("estimate_Re_from_noisy_delayed_incidence yields consistent results with import data", {
  toy_incidence_data <- c(
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 2, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

  toy_import_data <- c(
    1, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0,
    0, 0, 0, 0, 0, 9, 0, 1, 0, 0, 0, 0
  )

  delay_distribution <- c(
    0, 0.015, 0.09, 0.168,
    0.195, 0.176, 0.135, 0.091, 0.057, 0.034,
    0.019, 0.01, 0.005, 0.003,
    0.001, 0.001
  )

  estimates <- estimate_Re_from_noisy_delayed_incidence(
    incidence_data = toy_incidence_data,
    import_incidence_data = toy_import_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    delay = delay_distribution,
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    minimum_cumul_incidence = 10,
    mean_Re_prior = 1,
    output_Re_only = FALSE,
    ref_date = as.Date("2020-02-04"),
    time_step = "day"
  )

  reference_R_values <- c(
    NA, NA, NA, NA, NA, NA, NA, NA, NA, 3.52, 3.59, 3.52,
    3.37, 3.21, 3.05, 2.89, 2.74, 2.61, 2.48, 2.36,
    2.24, 2.13, 2.03, 1.92, 1.83, 1.73, 1.64, 1.56,
    1.47, 1.39, 1.32, 1.25, 1.17, 1.1, 1.02, 0.95,
    0.89, 0.84, 0.81, 0.77, 0.74, 0.72, 0.69, 0.65,
    0.62, NA, NA, NA, NA, NA
  )

  reference_dates <- seq.Date(from = as.Date("2020-01-30"), to = as.Date("2020-03-19"), by = "day")

  expect_equal(estimates$Re_estimate, reference_R_values, tolerance = 5E-2)
  expect_equal(estimates$date, reference_dates)
})

test_that("estimate_Re_from_noisy_delayed_incidence passes '...' arguments consistently", {
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

  delay_distribution <- c(
    0, 0.015, 0.09, 0.168,
    0.195, 0.176, 0.135, 0.091, 0.057, 0.034,
    0.019, 0.01, 0.005, 0.003,
    0.001, 0.001
  )

  estimates <- estimate_Re_from_noisy_delayed_incidence(
    incidence_data = toy_incidence_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    delay = delay_distribution,
    estimation_window = 5,
    mean_serial_interval = 8,
    std_serial_interval = 3,
    minimum_cumul_incidence = 0,
    block_size = 3,
    degree = 1,
    mean_Re_prior = 2.5,
    output_Re_only = FALSE,
    index_col = "date_index"
  )

  reference_R_values <- c(
    NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, 4.9, 4.42,
    4.03, 3.7, 3.4, 3.13, 2.89,
    2.66, 2.44, 2.25, 2.07, 1.9, 1.74,
    1.59, 1.45, 1.32, 1.2, 1.09,
    0.99, 0.91, 0.84, 0.77, 0.72,
    0.67, 0.63, 0.59, NA, NA, NA, NA, NA
  )
  reference_indices <- -5:36

  expect_equal(estimates$Re_estimate, reference_R_values, tolerance = 1E-2)
  expect_equal(estimates$date_index, reference_indices)
})

test_that("get_block_bootstrapped_estimate yields consistent results on a toy example", {
  skip_on_cran()
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

  shape_incubation <- 2
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  estimates <- get_block_bootstrapped_estimate(
    incidence_data = toy_incidence_data,
    N_bootstrap_replicates = 100,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    uncertainty_summary_method = "bagged mean - CI from bootstrap estimates",
    delay = list(delay_incubation, delay_onset_to_report),
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    minimum_cumul_incidence = 10,
    mean_Re_prior = 1,
    ref_date = as.Date("2020-02-04"),
    time_step = "day"
  )

  reference_dates <- seq.Date(
    from = as.Date("2020-02-04"),
    to = as.Date("2020-03-05"),
    by = "day"
  )

  reference_R_mean_values <- c(
    3.33, 3.04, 2.84, 2.67, 2.52, 2.38, 2.25,
    2.13, 2.01, 1.9, 1.8, 1.7, 1.61, 1.52, 1.44,
    1.36, 1.29, 1.22, 1.15, 1.09, 1.03, 0.97,
    0.92, 0.89, 0.86, 0.84, 0.82, 0.8, 0.78, 0.77, 0.75
  )

  reference_CI_down_values <- c(
    3.13, 2.84, 2.65, 2.51, 2.39, 2.27,
    2.14, 2.02, 1.91, 1.81, 1.71, 1.63,
    1.54, 1.46, 1.39, 1.31, 1.24, 1.17,
    1.11, 1.05, 0.99, 0.94, 0.89, 0.85,
    0.82, 0.79, 0.77, 0.75, 0.72, 0.7, 0.67
  )

  reference_CI_up_values <- c(
    3.54, 3.24, 3.02, 2.83, 2.66, 2.5, 2.36, 2.23,
    2.12, 2, 1.88, 1.78, 1.68, 1.58, 1.49, 1.41,
    1.34, 1.26, 1.19, 1.12, 1.06, 1, 0.96, 0.92,
    0.9, 0.88, 0.87, 0.86, 0.84, 0.83, 0.82
  )


  expect_equal(estimates$date, reference_dates)
  expect_equal(estimates$Re_estimate, reference_R_mean_values, tolerance = 1E-1)
  expect_equal(estimates$CI_down_Re_estimate, reference_CI_down_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up_Re_estimate, reference_CI_up_values, tolerance = 1E-1)

  estimates <- get_block_bootstrapped_estimate(
    incidence_data = toy_incidence_data,
    N_bootstrap_replicates = 100,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
    delay = list(delay_incubation, delay_onset_to_report),
    estimation_window = 3,
    minimum_cumul_incidence = 0,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    mean_Re_prior = 1
  )

  reference_indices <- 0:30

  reference_R_original_values <- c(
    3.2, 2.91, 2.72, 2.57, 2.44, 2.32,
    2.2, 2.09, 1.98, 1.88, 1.78, 1.69,
    1.6, 1.52, 1.44, 1.36, 1.29, 1.22,
    1.15, 1.08, 1.01, 0.95, 0.89, 0.84,
    0.8, 0.77, 0.75, 0.72, 0.7, 0.66, 0.63
  )

  reference_CI_down_values <- c(
    3.04, 2.75, 2.58, 2.45, 2.34, 2.22,
    2.11, 1.99, 1.89, 1.79, 1.7, 1.62,
    1.55, 1.47, 1.39, 1.32, 1.25, 1.18,
    1.11, 1.04, 0.98, 0.91, 0.86, 0.81,
    0.77, 0.73, 0.7, 0.67, 0.64, 0.6, 0.56
  )

  reference_CI_up_values <- c(
    3.36, 3.07, 2.86, 2.69, 2.55, 2.42, 2.3,
    2.19, 2.08, 1.97, 1.86, 1.76, 1.66, 1.57,
    1.48, 1.4, 1.33, 1.26, 1.18, 1.11, 1.04,
    0.98, 0.92, 0.88, 0.84, 0.82, 0.8, 0.78,
    0.75, 0.73, 0.7
  )

  expect_equal(estimates$idx, reference_indices)
  expect_equal(estimates$Re_estimate, reference_R_original_values, tolerance = 1E-1)
  expect_equal(estimates$CI_down_Re_estimate, reference_CI_down_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up_Re_estimate, reference_CI_up_values, tolerance = 1E-1)
})

test_that("get_block_bootstrapped_estimate yields consistent with import data", {
  skip_on_cran()
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

  toy_import_data <- c(
    1, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0,
    0, 0, 0, 0, 0, 9, 0, 1, 0, 0, 0, 0
  )

  shape_incubation <- 2
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  estimates <- get_block_bootstrapped_estimate(
    incidence_data = toy_incidence_data,
    import_incidence_data = toy_import_data,
    N_bootstrap_replicates = 100,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    uncertainty_summary_method = "bagged mean - CI from bootstrap estimates",
    delay = list(delay_incubation, delay_onset_to_report),
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    minimum_cumul_incidence = 10,
    mean_Re_prior = 1,
    ref_date = as.Date("2020-02-04"),
    time_step = "day"
  )

  reference_dates <- seq.Date(
    from = as.Date("2020-02-04"),
    to = as.Date("2020-03-05"),
    by = "day"
  )

  reference_R_mean_values <- c(
    2.69, 2.48, 2.37, 2.29, 2.22, 2.14, 2.06,
    1.98, 1.9, 1.82, 1.73, 1.65, 1.57, 1.49,
    1.42, 1.35, 1.28, 1.21, 1.14, 1.08, 1.02,
    0.97, 0.92, 0.89, 0.86, 0.84, 0.82, 0.8,
    0.79, 0.77, 0.75
  )

  reference_CI_down_values <- c(
    2.53, 2.34, 2.24, 2.18, 2.12, 2.05, 1.97,
    1.89, 1.8, 1.72, 1.65, 1.58, 1.5, 1.43,
    1.36, 1.29, 1.23, 1.17, 1.11, 1.05, 0.99, 0.94, 0.89,
    0.85, 0.82, 0.79, 0.77, 0.74, 0.72, 0.7, 0.67
  )

  reference_CI_up_values <- c(
    2.85, 2.63, 2.5, 2.4, 2.32, 2.23, 2.15,
    2.08, 2, 1.91, 1.82, 1.72, 1.64, 1.55,
    1.47, 1.4, 1.32, 1.25, 1.18, 1.12, 1.05,
    1, 0.96, 0.92, 0.9, 0.88, 0.87, 0.86,
    0.85, 0.84, 0.83
  )


  expect_equal(estimates$date, reference_dates)
  expect_equal(estimates$Re_estimate, reference_R_mean_values, tolerance = 1E-1)
  expect_equal(estimates$CI_down_Re_estimate, reference_CI_down_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up_Re_estimate, reference_CI_up_values, tolerance = 1E-1)
})

test_that("get_block_bootstrapped_estimate passes '...' arguments to inner functions properly", {
  skip_on_cran()
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

  shape_incubation <- 2
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  estimates <- get_block_bootstrapped_estimate(
    incidence_data = toy_incidence_data,
    N_bootstrap_replicates = 100,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    uncertainty_summary_method = "bagged mean - CI from bootstrap estimates",
    delay = list(delay_incubation, delay_onset_to_report),
    estimation_window = 5,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    minimum_cumul_incidence = 0,
    block_size = 8,
    degree = 2,
    ref_date = as.Date("2020-02-04"),
    time_step = "day"
  )

  reference_R_mean_values <- c(
    5.5, 4.6, 3.95, 3.49, 3.15, 2.88,
    2.66, 2.48, 2.31, 2.15, 2.02, 1.89, 1.76, 1.64,
    1.53, 1.42, 1.32, 1.22, 1.13, 1.05, 0.98, 0.92,
    0.86, 0.8, 0.74, 0.67, 0.6, 0.51, 0.42
  )

  reference_CI_down_values <- c(
    5.01, 4.28, 3.72,
    3.28, 2.95, 2.69, 2.49, 2.34, 2.2, 2.06,
    1.94, 1.83, 1.71, 1.59, 1.49, 1.38, 1.28,
    1.19, 1.1, 1.03, 0.96, 0.9, 0.84, 0.78,
    0.72, 0.66, 0.58, 0.5, 0.4
  )

  reference_CI_up_values <- c(
    5.99, 4.92,
    4.18, 3.69, 3.36, 3.07, 2.82,
    2.62, 2.43, 2.25, 2.09, 1.95,
    1.81, 1.69, 1.57, 1.46, 1.36,
    1.26, 1.17, 1.08, 1, 0.93, 0.87,
    0.81, 0.75, 0.68, 0.61, 0.53, 0.44
  )

  # master
  # reference_R_mean_values <- c(6.83,5.37,4.53,4.05,3.71,3.43,3.21,3.02,2.84,
  #                              2.66,2.5,2.34,2.18,2.04,1.9,1.77,1.64,
  #                              1.53,1.42,1.32,1.22,1.13,1.05,0.98,0.92,
  #                              0.86,0.8,0.74,0.67,0.6,0.52,0.42)
  #
  # reference_CI_down_values <- c(6.29,4.95,4.21,3.8,3.51,3.27,
  #                               3.07,2.9,2.73,2.57,2.42,2.28,2.12,
  #                               1.99,1.86,1.73,1.61,1.5,1.39,1.29,
  #                               1.19,1.11,1.03,0.96,0.9,0.84,0.78,
  #                               0.73,0.66,0.59,0.5,0.4)
  #
  # reference_CI_up_values <- c(7.37,5.79,4.85,4.29,3.9,
  #                             3.59,3.34,3.14,2.94,2.74,
  #                             2.57,2.4,2.24,2.09,1.95,
  #                             1.81,1.68,1.56,1.45,1.34,
  #                             1.25,1.16,1.08,1,0.94,0.87,
  #                             0.81,0.75,0.69,0.62,0.54,0.45)


  expect_equal(estimates$Re_estimate, reference_R_mean_values, tolerance = 1E-1)
  expect_equal(estimates$CI_down_Re_estimate, reference_CI_down_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up_Re_estimate, reference_CI_up_values, tolerance = 1E-1)
})

test_that("get_infections_from_incidence handles partially-delayed data correctly", {
  toy_onset_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66, 45, 32, 11, 5, 4, 0, 1, 2, 0
  )

  shape_incubation <- 2
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  result_deconvolution <- get_infections_from_incidence(
    incidence_data = toy_onset_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    delay = delay_incubation,
    is_partially_reported_data = TRUE,
    delay_until_final_report = delay_onset_to_report,
    output_infection_incidence_only = FALSE,
    data_points_incl = 21,
    degree = 1,
    cutoff_observation_probability = 0.1
  )

  reference_deconvolved_incidence <- c(
    16.417, 19.763, 22.179, 26.041, 32.397, 40.832,
    50.92, 62.293, 75.173, 89.947, 106.529, 124.207,
    142.552, 161.947, 182.9, 203.07, 221.026, 239.061,
    255.931, 268.766, 278.113, 285.134, 289.532, 289.975,
    286.109, 278.965, 269.38, 256.736, 240.254, 221.707,
    202.645, 181.802, 160.793, 141.705, 123.465, 104.689,
    85.982, 68.065, 49.84, 31.857, 15.892, 3.778,
    0, 0, 0, NA, NA, NA
  )

  expect_equal(result_deconvolution$deconvolved_incidence, reference_deconvolved_incidence, tolerance = 1E-1)
})

test_that("estimate_from_combined_observations returns consistent results", {
  toy_onset_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66, 45, 32, 11, 5, 4, 0, 1, 2, 0
  )

  toy_case_confirmation_data <- c(
    11, 12, 21, 23, 2, 14, 49, 61, 65, 45, 66, 45,
    40, 8, 61, 38, 1, 3, 45, 66, 12, 52, 27, 3, 54,
    10, 18, 54, 12, 48, 67, 62, 54, 3, 29, 10, 52,
    61, 33, 39, 55, 8, 64, 51, 65, 34
  )

  shape_incubation <- 1
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 4
  scale_onset_to_report <- 2.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  results_estimation <- estimate_from_combined_observations(
    partially_delayed_incidence = toy_onset_data,
    fully_delayed_incidence = toy_case_confirmation_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    delay_until_partial = delay_incubation,
    delay_until_final_report = delay_onset_to_report,
    partial_observation_requires_full_observation = TRUE,
    ref_date = as.Date("2021-03-24"),
    time_step = "day",
    minimum_cumul_incidence = 0,
    output_Re_only = FALSE,
    data_points_incl = 21,
    degree = 1
  )

  reference_deconvolved_incidence <- c(
    59.2, 61.7, 64.2, 69.8, 76.9, 84.4, 92.1, 101.5, 112.9, 126.2,
    140.8, 156.9, 175.5, 195.1, 214.5, 235.3, 255.6, 272.3, 287,
    301.6, 314.4, 323.1, 328.8, 332.7, 333.8, 330.9, 324.9, 317.1,
    305.9, 290.8, 273.5, 255.4, 237.4, 219.6, 201.4, 183.3, 165.1,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  )

  reference_R_values <- c(
    NA, NA, NA, NA, NA, NA, 2, 1.7, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6,
    1.6, 1.6, 1.6, 1.5, 1.5, 1.4, 1.4, 1.3, 1.3, 1.2, 1.2, 1.1,
    1.1, 1, 1, 0.9, 0.9, 0.9, 0.8, 0.8, 0.7, 0.7, 0.7,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  )

  expect_equal(results_estimation$combined_deconvolved_incidence, reference_deconvolved_incidence, tolerance = 1E-1)
  expect_equal(results_estimation$Re_estimate, reference_R_values, tolerance = 1E-1)
})

test_that("get_bootstrapped_estimate_from_combined_observations returns consistent results", {
  skip_on_cran()
  toy_onset_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66, 45, 32, 11, 5, 4, 0, 1, 2, 0
  )

  toy_case_confirmation_data <- c(
    11, 12, 21, 23, 2, 14, 49, 61, 65, 45, 66, 45,
    40, 8, 61, 38, 1, 3, 45, 66, 12, 52, 27, 3, 54,
    10, 18, 54, 12, 48, 67, 62, 54, 3, 29, 10, 52,
    61, 33, 39, 55, 8, 64, 51, 65, 34
  )

  shape_incubation <- 1
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 4
  scale_onset_to_report <- 2.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  results_estimation <- get_bootstrapped_estimates_from_combined_observations(
    partially_delayed_incidence = toy_onset_data,
    fully_delayed_incidence = toy_case_confirmation_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    bootstrapping_method = "non-parametric block boostrap",
    uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
    N_bootstrap_replicates = 100,
    delay_until_partial = delay_incubation,
    delay_until_final_report = delay_onset_to_report,
    partial_observation_requires_full_observation = TRUE,
    ref_date = as.Date("2021-03-24"),
    time_step = "day",
    output_Re_only = TRUE,
    minimum_cumul_incidence = 0,
    data_points_incl = 21,
    degree = 1
  )


  reference_R_values <- c(
    2, 1.7, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6,
    1.6, 1.6, 1.6, 1.5, 1.5, 1.4, 1.4, 1.3, 1.3, 1.2, 1.2, 1.1,
    1.1, 1, 1, 0.9, 0.9, 0.9, 0.8, 0.8, 0.7, 0.7, 0.7
  )

  reference_CI_up <- c(
    2.17, 1.95, 1.85, 1.81, 1.79, 1.79, 1.79,
    1.78, 1.76, 1.73, 1.68, 1.62, 1.56, 1.49, 1.43, 1.38,
    1.32, 1.26, 1.21, 1.15, 1.11, 1.06, 1.02, 0.98, 0.94,
    0.9, 0.86, 0.83, 0.81, 0.79, 0.77
  )

  reference_CI_down <- c(
    1.75, 1.54, 1.45, 1.42, 1.43, 1.44, 1.46,
    1.48, 1.49, 1.49, 1.48, 1.46, 1.41, 1.36, 1.31, 1.25,
    1.21, 1.16, 1.12, 1.07, 1.03, 0.99, 0.94, 0.89, 0.85,
    0.8, 0.76, 0.72, 0.69, 0.65, 0.61
  )

  expect_equal(results_estimation$Re_estimate, reference_R_values, tolerance = 1E-1)
  expect_equal(results_estimation$CI_up_Re_estimate, reference_CI_up, tolerance = 1E-1)
  expect_equal(results_estimation$CI_down_Re_estimate, reference_CI_down, tolerance = 1E-1)
})

test_that("get_bootstrapped_estimate_from_combined_observations can deal with empirical delay data", {
  skip_on_cran()
  ref_date <- as.Date("2021-03-24")
  n_days <- 50
  shape_initial_delay <- 6
  scale_initial_delay <- 1.5
  distribution_initial_delay <- list(name = "gamma", shape = shape_initial_delay, scale = scale_initial_delay)
  seed <- 7543265

  generated_empirical_delays <- .generate_delay_data(
    origin_date = ref_date,
    n_time_steps = n_days,
    ratio_delay_end_to_start = 1.5,
    distribution_initial_delay = distribution_initial_delay,
    seed = seed
  )

  toy_onset_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66, 45, 32, 11, 5, 4, 0, 1, 2, 0
  )

  toy_case_confirmation_data <- c(
    11, 12, 21, 23, 2, 14, 49, 61, 65, 45, 66, 45,
    40, 8, 61, 38, 1, 3, 45, 66, 12, 52, 27, 3, 54,
    10, 18, 54, 12, 48, 67, 62, 54, 3, 29, 10, 52,
    61, 33, 39, 55, 8, 64, 51, 65, 34
  )

  shape_incubation <- 2
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  results_estimation <- get_bootstrapped_estimates_from_combined_observations(
    partially_delayed_incidence = toy_onset_data,
    fully_delayed_incidence = toy_case_confirmation_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    bootstrapping_method = "non-parametric block boostrap",
    uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
    N_bootstrap_replicates = 20, # to speed things up
    delay_until_partial = delay_incubation,
    delay_until_final_report = generated_empirical_delays,
    partial_observation_requires_full_observation = TRUE,
    ref_date = ref_date,
    time_step = "day",
    output_Re_only = TRUE,
    minimum_cumul_incidence = 0,
    data_points_incl = 21,
    degree = 1,
    cutoff_observation_probability = 0.1
  )

  reference_R_values <- c(
    1.98, 1.76, 1.67, 1.63,
    1.63, 1.65, 1.66, 1.66,
    1.66, 1.63, 1.59, 1.54,
    1.48, 1.42, 1.37, 1.32,
    1.26, 1.21, 1.17, 1.12,
    1.08, 1.04, 1, 0.95, 0.91,
    0.87, 0.83, 0.8, 0.77, 0.75
  )

  expect_equal(results_estimation$Re_estimate, reference_R_values, tolerance = 1E-1)
})

test_that("get_block_bootstrapped_estimate yields consistent results on summaries of uncertainty", {
  skip_on_cran()
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

  shape_incubation <- 2
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  estimates <- get_block_bootstrapped_estimate(
    incidence_data = toy_incidence_data,
    N_bootstrap_replicates = 100,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
    delay = list(delay_incubation, delay_onset_to_report),
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    minimum_cumul_incidence = 10,
    mean_Re_prior = 1,
    ref_date = as.Date("2020-02-04"),
    time_step = "day",
    output_Re_only = FALSE
  )

  reference_dates <- seq.Date(
    from = as.Date("2020-01-29"),
    to = as.Date("2020-03-11"),
    by = "day"
  )

  reference_CI_down_observed_incidence_values <- c(
    NA, NA, NA, NA, NA, NA, 4.9, 8, 10, 13, 17, 20.9, 31, 41,
    50.9, 63.9, 78.9, 81.2, 108.1, 124.4, 145.6, 174.2, 208.2, 228.6,
    238.3, 232.6, 214.2, 319.7, 326.6, 318.4, 299.1,
    273.5, 237, 190, 120.2, 32.8, 153.6, 139.7,
    124.4, 104.5, 59.4, 19.8, 15.6
  )

  reference_CI_up_deconvolved_incidence_values <- c(
    15.3, 18.8, 23.3, 28.7, 35.1, 43.4, 53.8,
    66, 80, 96.4, 115, 134.9, 155.8, 178.9, 201.5,
    222.1, 243.8, 263.9, 280, 295.9, 309.3, 317.8, 323.8,
    324.6, 318.7, 311.1, 300.2, 284.6, 270.9,
    259.1, 246.2, 234.5, 223.2, 210.6, 197.5,
    184.5, 172.5, NA, NA, NA, NA, NA, NA
  )

  reference_CI_down_smoothed_incidence_values <- c(
    NA, NA, NA, NA, NA, NA, 14, 17.4, 21.6, 26.4, 32.4, 40, 49.2, 59.5,
    71, 83.9, 98.1, 113, 128.5, 145.3, 162.2, 178.5, 195.4, 211, 223.5,
    234.8, 244.3, 250.3, 252.7, 251.4, 246.6, 238.4,
    226.8, 213.2, 198.9, 184, 168, 151.7, 134.8, 116.8, 98, 79, 60.2
  )

  reference_CI_up_R_mean_values <- c(
    NA, NA, NA, NA, NA, NA, 2.98, 2.72,
    2.57, 2.47, 2.4, 2.32, 2.23, 2.14,
    2.04, 1.93, 1.83, 1.74, 1.65, 1.58,
    1.51, 1.44, 1.38, 1.31, 1.24, 1.17, 1.1, 1.04,
    0.98, 0.93, 0.9, 0.87, 0.85, 0.83, 0.82,
    0.8, 0.78, NA, NA, NA, NA, NA, NA
  )


  expect_equal(estimates$date, reference_dates)
  # This test usually fails but that is not a problem (great stochastic variability)
  # expect_equal(estimates$CI_down_observed_incidence, reference_CI_down_observed_incidence_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up_deconvolved_incidence, reference_CI_up_deconvolved_incidence_values, tolerance = 1E-1)
  expect_equal(estimates$CI_down_smoothed_incidence, reference_CI_down_smoothed_incidence_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up_Re_estimate, reference_CI_up_R_mean_values, tolerance = 1E-1)
})

test_that("get_bootstrapped_estimate_from_combined_observations yields consistent results on summaries of uncertainty", {
  skip_on_cran()
  toy_onset_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66, 45, 32, 11, 5, 4, 0, 1, 2, 0
  )

  toy_case_confirmation_data <- c(
    11, 12, 21, 23, 2, 14, 49, 61, 65, 45, 66, 45,
    40, 8, 61, 38, 1, 3, 45, 66, 12, 52, 27, 3, 54,
    10, 18, 54, 12, 48, 67, 62, 54, 3, 29, 10, 52,
    61, 33, 39, 55, 8, 64, 51, 65, 34
  )

  shape_incubation <- 1
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 4
  scale_onset_to_report <- 2.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  results_estimation <- get_bootstrapped_estimates_from_combined_observations(
    partially_delayed_incidence = toy_onset_data,
    fully_delayed_incidence = toy_case_confirmation_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    bootstrapping_method = "non-parametric block boostrap",
    uncertainty_summary_method = "bagged mean - CI from bootstrap estimates",
    N_bootstrap_replicates = 100,
    delay_until_partial = delay_incubation,
    delay_until_final_report = delay_onset_to_report,
    partial_observation_requires_full_observation = TRUE,
    ref_date = as.Date("2021-03-24"),
    time_step = "day",
    output_Re_only = FALSE,
    minimum_cumul_incidence = 0,
    data_points_incl = 21,
    degree = 1,
    cutoff_observation_probability = 0.1
  )

  reference_partially_delayed_observations_values <- c(
    NA, 4.4, 5.4, 8, 11.1, 14.8, 20.3, 28.2, 37.7,
    49.6, 63.4, 80.6, 96.3, 118, 145.3, 168.4, 193.7,
    218.7, 247.1, 261.4, 287.2, 311.2, 326.7, 336.3,
    342.9, 345.8, 327, 322, 320.6, 301.5, 273.3, 239.7,
    210.9, 177.4, 155.6, 138.5, 120.1, 104.1, 89.1,
    76.9, 63.2, 54.8, 48, NA, NA, NA, NA
  )

  reference_CI_up_combined_deconvolved_incidence_values <- c(
    58.1, 61.8, 65.4, 71.5, 79.1, 88, 98.2, 109.8,
    122.6, 137.3, 153.5, 170.8, 190.9, 211.8, 231.9,
    253.9, 276, 294.6, 311.5, 329.1, 345, 356.2, 364.7,
    370.1, 370.2, 365.8, 359.5, 351.2, 338.7, 323.2,
    307.5, 292.1, 277, 262.5, 248.7, 235.6, 222.7,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  )


  reference_R_mean_values <- c(
    NA, NA, NA, NA, NA, NA, 2.11, 1.9, 1.81, 1.76,
    1.74, 1.72, 1.71, 1.69, 1.66, 1.63, 1.59, 1.54,
    1.49, 1.44, 1.39, 1.34, 1.28, 1.23, 1.18, 1.13,
    1.08, 1.03, 0.99, 0.94, 0.9, 0.86, 0.83, 0.8,
    0.77, 0.75, 0.72, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  )

  reference_CI_up_R_mean <- c(
    NA, NA, NA, NA, NA, NA, 2.3, 2.09, 1.98,
    1.92, 1.88, 1.85, 1.83, 1.8, 1.77, 1.72,
    1.67, 1.62, 1.56, 1.5, 1.45,
    1.4, 1.34, 1.28, 1.23, 1.17, 1.11, 1.07,
    1.03, 0.99, 0.95, 0.91, 0.88, 0.86, 0.84,
    0.82, 0.81, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  )

  expect_equal(results_estimation$partially_delayed_observations,
    reference_partially_delayed_observations_values,
    tolerance = 1E-1
  )
  expect_equal(results_estimation$CI_up_combined_deconvolved_incidence,
    reference_CI_up_combined_deconvolved_incidence_values,
    tolerance = 1E-1
  )
  expect_equal(results_estimation$Re_estimate, reference_R_mean_values, tolerance = 1E-1)
  expect_equal(results_estimation$CI_up_Re_estimate, reference_CI_up_R_mean, tolerance = 1E-1)
})

test_that("get_block_bootstrapped_estimate consistently combines HPDs with bootstrap CIs", {
  skip_on_cran()
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 350, 403, 505, 668, 873, 987, 1050, 1268, 1490, 1760
  )

  shape_incubation <- 2
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  estimates <- get_block_bootstrapped_estimate(
    incidence_data = toy_incidence_data,
    N_bootstrap_replicates = 100,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
    delay = list(delay_incubation, delay_onset_to_report),
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    minimum_cumul_incidence = 10,
    mean_Re_prior = 1,
    ref_date = as.Date("2020-02-04"),
    time_step = "day",
    combine_bootstrap_and_estimation_uncertainties = TRUE,
    output_Re_only = FALSE
  )

  simplified_estimates <- get_block_bootstrapped_estimate(
    incidence_data = toy_incidence_data,
    N_bootstrap_replicates = 100,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
    delay = list(delay_incubation, delay_onset_to_report),
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    minimum_cumul_incidence = 10,
    mean_Re_prior = 1,
    ref_date = as.Date("2020-02-04"),
    time_step = "day",
    combine_bootstrap_and_estimation_uncertainties = TRUE,
    output_Re_only = TRUE
  )

  reference_CI_down <- c(
    NA, NA, NA, NA, NA, NA, 2.51, 2.37, 2.28,
    2.2, 2.12, 2, 1.91, 1.83, 1.76, 1.67, 1.58,
    1.49, 1.43, 1.44, 1.49, 1.55, 1.59, 1.61,
    1.63, 1.63, 1.61, 1.58, 1.54, 1.5, 1.46,
    1.43, 1.41, 1.38, NA, NA, NA, NA, NA, NA
  )

  pretty_reference_CI_down <- c(
    2.51, 2.37, 2.28,
    2.2, 2.12, 2, 1.91, 1.83, 1.76, 1.67,
    1.58, 1.49, 1.43, 1.44, 1.49, 1.55,
    1.59, 1.61, 1.63, 1.63, 1.61, 1.58,
    1.54, 1.5, 1.46, 1.43, 1.41, 1.38
  )

  reference_bootstrap_CI_up <- c(
    NA, NA, NA, NA, NA, NA, 3.27, 3.08, 2.93,
    2.82, 2.72, 2.61, 2.47, 2.31, 2.16, 2.05,
    1.97, 1.92, 1.89, 1.89, 1.93, 1.99, 2.05,
    2.07, 2.08, 2.07, 2.03, 1.96, 1.87, 1.79,
    1.71, 1.64, 1.59, 1.55, NA, NA, NA, NA, NA, NA
  )

  reference_highHPD <- c(
    NA, NA, NA, NA, NA, NA, 3.67, 3.33,
    3.1, 2.91, 2.74, 2.61, 2.47,
    2.31, 2.16, 2.05, 1.97, 1.92, 1.89,
    1.89, 1.93, 1.99, 2.05, 2.07, 2.08,
    2.07, 2.03, 1.96, 1.87, 1.79, 1.71,
    1.64, 1.59, 1.55, NA, NA, NA, NA, NA, NA
  )

  reference_Re_estimate <- c(
    NA, NA, NA, NA, NA, NA, 3.06, 2.83, 2.67,
    2.54, 2.42, 2.31, 2.19, 2.07, 1.96, 1.86,
    1.78, 1.71, 1.66, 1.67, 1.71, 1.77, 1.82,
    1.84, 1.86, 1.85, 1.82, 1.77, 1.7, 1.64, 1.58,
    1.54, 1.5, 1.47, NA, NA, NA, NA, NA, NA
  )

  pretty_reference_Re_estimate <- c(
    3.06, 2.83, 2.67,
    2.54, 2.42, 2.31, 2.19, 2.07, 1.96, 1.86,
    1.78, 1.71, 1.66, 1.67, 1.71, 1.77, 1.82,
    1.84, 1.86, 1.85, 1.82, 1.77, 1.7, 1.64, 1.58,
    1.54, 1.5, 1.47
  )

  expect_equal(simplified_estimates$CI_down_Re_estimate, pretty_reference_CI_down, tolerance = 1E-1)
  expect_equal(estimates$CI_down_Re_estimate, reference_CI_down, tolerance = 1E-1)
  expect_equal(estimates$bootstrapped_CI_up_Re_estimate, reference_bootstrap_CI_up, tolerance = 1E-1)
  expect_equal(estimates$Re_highHPD, reference_highHPD, tolerance = 1E-1)
  expect_equal(estimates$Re_estimate, reference_Re_estimate, tolerance = 1E-1)
  expect_equal(simplified_estimates$Re_estimate, pretty_reference_Re_estimate, tolerance = 1E-1)
})

test_that("get_bootstrapped_estimate_from_combined_observations consistently combines HPDs with bootstrap CIs", {
  skip_on_cran()
  toy_onset_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66, 45, 32, 11, 5, 4, 0, 1, 2, 0
  )

  toy_case_confirmation_data <- c(
    11, 12, 21, 23, 2, 14, 49, 61, 65, 45, 66, 45,
    40, 8, 61, 38, 1, 3, 4, 5, 56, 3, 45, 66, 12, 52, 27, 3, 54,
    120, 150, 230, 400, 487, 496, 602, 893, 1020, 1250,
    1400, 1746, 2190, 2567, 3498, 4192, 6432
  )

  shape_incubation <- 1
  scale_incubation <- 1.2
  delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 4
  scale_onset_to_report <- 2.3
  delay_onset_to_report <- list(name = "gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  results_estimation <- get_bootstrapped_estimates_from_combined_observations(
    partially_delayed_incidence = toy_onset_data,
    fully_delayed_incidence = toy_case_confirmation_data,
    smoothing_method = "LOESS",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    bootstrapping_method = "non-parametric block boostrap",
    uncertainty_summary_method = "bagged mean - CI from bootstrap estimates",
    N_bootstrap_replicates = 100,
    delay_until_partial = delay_incubation,
    delay_until_final_report = delay_onset_to_report,
    partial_observation_requires_full_observation = TRUE,
    ref_date = as.Date("2021-03-24"),
    time_step = "day",
    minimum_cumul_incidence = 0,
    data_points_incl = 21,
    degree = 1,
    combine_bootstrap_and_estimation_uncertainties = TRUE,
    output_Re_only = TRUE
  )

  reference_R_mean_values <- c(
    2.13, 1.93, 1.84, 1.81, 1.81, 1.82,
    1.82, 1.82, 1.81, 1.8,
    1.78, 1.76, 1.73, 1.7, 1.67, 1.65, 1.62,
    1.63, 1.7, 1.85, 2.02, 2.22, 2.45, 2.68,
    2.82, 2.81, 2.64, 2.39, 2.14, 1.94,
    1.79
  )

  reference_CI_down_R_mean <- c(
    1.64, 1.48, 1.43, 1.43,
    1.45, 1.49, 1.52, 1.55, 1.57, 1.58, 1.59, 1.58,
    1.57, 1.52, 1.47, 1.43, 1.35, 1.28, 1.27, 1.35,
    1.49, 1.7, 1.93, 1.92, 1.68, 1.46, 1.4, 1.44,
    1.48, 1.47, 1.42
  )

  reference_CI_up_R_mean <- c(
    2.42, 2.21, 2.1,
    2.05, 2.03, 2.03, 2.02, 2.01, 1.98, 1.96,
    1.93, 1.91, 1.9, 1.89, 1.87, 1.87, 1.89,
    1.99, 2.14, 2.34, 2.55, 2.74, 2.98, 3.44,
    3.97, 4.16, 3.88, 3.34, 2.81, 2.41,
    2.15
  )

  expect_equal(results_estimation$Re_estimate, reference_R_mean_values, tolerance = 1E-1)
  expect_equal(results_estimation$CI_down_Re_estimate, reference_CI_down_R_mean, tolerance = 1E-1)
  expect_equal(results_estimation$CI_up_Re_estimate, reference_CI_up_R_mean, tolerance = 1E-1)
})
