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

  reference_values <- c(4.9,6.6,8.3,10.9,14.4,18.8,27,36.5,
                        47.4,60.3,75.1,91.9,110.9,133,157.6,
                        185.9,216.9,246.6,272.5,299.2,320.7,
                        333.6,343,345.4,336.7,323.2,306.2,283,
                        261.1,239.6,216.6,193.9,171.2,145.3,
                        117.2,89.1,63.4)
  reference_offset <- -5

  expect_lte(max(abs(.get_values(deconvolved_incidence) - reference_values)), expected = 1)
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
    delay = list(
      distribution_incubation,
      generated_empirical_delays
    ),
    ref_date = ref_date,
    min_number_cases = 20,
    threshold_chi_squared = 1,
    max_iterations = 100,
    fit = "none"
  )

  reference_values <- c(3.6,4.9,6.2,8.2,11,14.8,
                        21.8,30.2,40.4,53.2,69,87.6,
                        108.8,133.8,162.3,195.9,234.4,
                        273.5,309,344.4,372.2,
                        388,396.9,394.8,377.6,353.2,
                        323.9,287.9,254.6,
                        223.8,194.3,167.3,141.2,
                        113.8,87.6,64.1,43.9)

  reference_offset <- -11

  expect_lte(max(abs(.get_values(deconvolved_incidence) - reference_values)), expected = 0.1)
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
    delay = list(
      distribution_incubation,
      generated_empirical_delays
    ),
    ref_date = ref_date,
    threshold_chi_squared = 1,
    max_iterations = 100,
    fit = fit,
    min_number_cases = min_number_cases
  )

  reference_values <- c(2.4,3.3,4.3,5.7,7.8,10.5,15.6,21.9,29.7,39.8,52.6,
                        68.5,87.7,111.4,138.9,171.6,208.8,246.8,282.8,321,
                        354,376.6,393.3,399.8,390.7,373.6,350.6,319,288.1,
                        257.5,225.5,194.9,165.8,135.6,105.6,77.7,53.6)

  reference_offset <- -14

  expect_lte(max(abs(.get_values(deconvolved_incidence) - reference_values)), expected = 0.1) # Work-around because weird behaviour of expect_equal()
  expect_identical(.get_offset(deconvolved_incidence), reference_offset)
})

test_that("deconvolve_incidence handles incidence with 0s at the end", {
  toy_incidence_data <- c(
    6, 8, 10, 13, 17, 22, 31, 41, 52, 65, 80, 97, 116,
    138, 162, 189, 218, 245, 268, 292, 311, 322, 330,
    332, 324, 312, 297, 276, 256, 236, 214, 192, 170,
    145, 118, 91, 66, 32, 14, 0, 0
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

  reference_values <- c(4.8,6.4,8.1,10.7,14.1,18.4,
                        26.4,35.8,46.9,60.1,74.9,91.5,
                        110.1,132,156.4,184.6,216,246.4,
                        273,300.4,322.3,335.8,345.5,
                        347.9,339,325.1,307.4,283.6,
                        261.5,240,217.3,195.1,171.8,
                        144,112.1,79.1,48.9,18.2,5.2,0,0)
  reference_offset <- -5

  expect_lte(max(abs(.get_values(deconvolved_incidence) - reference_values)), expected = 1)
  expect_identical(.get_offset(deconvolved_incidence), reference_offset)
})
