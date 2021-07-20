test_that("deconvolve_incidence yields consistent results on a toy constant-delay example", {
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

  deconvolved_incidence <- deconvolve_incidence(
    incidence_data = toy_incidence_data,
    deconvolution_method = "Richardson-Lucy delay distribution",
    delay_incubation = delay_distribution,
    threshold_chi_squared = 1,
    max_iterations = 100
  )

  reference_values <- c(
    6, 8, 9, 12, 15, 19, 27, 37, 47, 60, 75, 92, 111, 133,
    158, 186, 217, 247, 272, 299, 321, 334, 343, 345, 337,
    323, 306, 283, 261, 239, 216, 193, 169, 143, 114, 85, 59
  )
  reference_offset <- -5

  expect_equal(.get_values(deconvolved_incidence), reference_values, tolerance = 1)
  expect_identical(.get_offset(deconvolved_incidence), reference_offset)
})

test_that("deconvolve_incidence yields consistent results on a toy moving-through-time example", {
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

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

  shape_incubation <- 3.2
  scale_incubation <- 1.3
  distribution_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

  deconvolved_incidence <- deconvolve_incidence(
    incidence_data = toy_incidence_data,
    deconvolution_method = "Richardson-Lucy delay distribution",
    delay_incubation = distribution_incubation,
    delay_onset_to_report = generated_empirical_delays,
    ref_date = ref_date,
    min_number_cases = 20,
    threshold_chi_squared = 1,
    max_iterations = 100
  )

  reference_values <- c(3.53,4.81,6.09,8.01,10.72,
                        14.37,21.12,29.24,39.02,
                        51.75,67.99,87.85,110.71,
                        136.73,164.73,196.98,234.54,
                        273.46,309.18,344.9,372.96,
                        389.07,398.1,395.98,378.43,
                        353.4,323.17,286.36,252.26,
                        220.97,191.15,163.87,137.36,
                        109.48,83.23,60.16,40.46)
  reference_offset <- -11

  expect_lte( max(abs(.get_values(deconvolved_incidence) - reference_values)), expected = 0.1)
  expect_identical(.get_offset(deconvolved_incidence), reference_offset)
})

test_that("deconvolve_incidence takes into account extra parameters for get_matrix_empirical_delay_distr", {
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66
  )

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
    ratio_delay_end_to_start = 3,
    distribution_initial_delay = distribution_initial_delay,
    seed = seed
  )

  shape_incubation <- 3.2
  scale_incubation <- 1.3
  distribution_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)
  fit <- "gamma"
  min_number_cases <- 20

  deconvolved_incidence <- deconvolve_incidence(
    incidence_data = toy_incidence_data,
    deconvolution_method = "Richardson-Lucy delay distribution",
    delay_incubation = distribution_incubation,
    delay_onset_to_report = generated_empirical_delays,
    ref_date = ref_date,
    threshold_chi_squared = 1,
    max_iterations = 100,
    fit = fit,
    min_number_cases = min_number_cases
  )

  reference_values <- c(2.44,3.32,4.26,5.72,7.77,
                        10.52,15.61,21.91,29.68,39.79,
                        52.62,68.5,87.73,111.35,138.94,
                        171.69,209.07,247.21,283.37,321.86,
                        355.1,377.96,394.83,401.25,391.94,
                        374.41,350.6,318.05,285.99,254.04,
                        220.78,188.94,158.86,128.11,98.12,
                        70.87,47.9)

  reference_offset <- -14

  expect_lte( max(abs(.get_values(deconvolved_incidence) - reference_values)), expected = 0.1 ) # Work-around because weird behaviour of expect_equal()
  expect_identical(.get_offset(deconvolved_incidence), reference_offset)
})
