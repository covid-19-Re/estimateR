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
    delay = delay_distribution,
    threshold_chi_squared = 1,
    max_iterations = 100
  )

  reference_values <- c(4.9,6.6,8.3,10.9,14.4,18.8,27,36.5,47.4,
                        60.3,75.1,91.9,110.9,133,157.6,185.9,216.9,
                        246.6,272.5,299.2,320.7,333.6,343,345.4,
                        336.7,323.2,306.1,282.8,260.8,239.2,215.9,
                        192.7,169.4,142.8,113.9,85.3,59.4)
  reference_offset <- -5

  expect_lte( max(abs(.get_values(deconvolved_incidence) - reference_values)), expected = 1)
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
    delay = list(distribution_incubation,
                  generated_empirical_delays),
    ref_date = ref_date,
    min_number_cases = 20,
    threshold_chi_squared = 1,
    max_iterations = 100,
    fit = "none"
  )

 reference_values <- c(3.6,4.9,6.2,8.2,11,14.8,
                       21.8,30.2,40.4,53.2,69,
                       87.6,108.7,133.8,162.3,
                       195.9,234.5,273.8,309.5,
                       345.2,373.3,389.3,398.3,
                       396.1,378.5,353.5,323.2,
                       286.4,252.3,221,191.1,
                       163.9,137.3,109.5,83.2,
                       60.2,40.5)

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
    delay = list(distribution_incubation,
                  generated_empirical_delays),
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
