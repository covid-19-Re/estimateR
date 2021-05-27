test_that("smooth_deconvolve_estimate yields consistent results on a toy example", {

  toy_incidence_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                          138,162,189,218,245,268,292,311,322,330,
                          332,324,312,297,276,256,236,214,192,170,
                          145,118,91,66)

  delay_distribution <- c(0,0.015,0.09,0.168,
                          0.195,0.176,0.135,0.091,0.057,0.034,
                          0.019,0.01,0.005,0.003,
                          0.001,0.001)

  estimates <- smooth_deconvolve_estimate(incidence_data = toy_incidence_data,
                                          smoothing_method = "LOESS",
                                          deconvolution_method = "Richardson-Lucy delay distribution",
                                          estimation_method = "EpiEstim sliding window",
                                          delay_incubation = delay_distribution,
                                          estimation_window = 3,
                                          mean_serial_interval = 4.8,
                                          std_serial_interval  = 2.3,
                                          minimum_cumul_incidence = 10,
                                          mean_Re_prior = 1,
                                          output_Re_only = FALSE,
                                          ref_date = as.Date("2020-02-04"),
                                          time_step = "day")

  reference_R_values <- c(NA,NA,NA,NA,4.785,3.426,2.882,2.644,2.518,2.426,2.343,
                          2.255,2.16,2.064,1.969,1.873,1.78,1.691,1.605,1.522,1.444,
                          1.37,1.3,1.23,1.159,1.09,1.023,0.96,0.903,0.855,0.818,0.788,
                          0.762,0.737,0.709,0.677,0.641,NA,NA,NA,NA,NA)

  reference_dates <- seq.Date(from = as.Date("2020-01-30"), to = as.Date("2020-03-11"), by = "day")

  expect_equal(estimates$R_mean, reference_R_values, tolerance = 1E-2)
  expect_equal(estimates$date, reference_dates)
})

test_that("smooth_deconvolve_estimate passes '...' arguments consistently", {



  toy_incidence_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                          138,162,189,218,245,268,292,311,322,330,
                          332,324,312,297,276,256,236,214,192,170,
                          145,118,91,66)

  delay_distribution <- c(0,0.015,0.09,0.168,
                          0.195,0.176,0.135,0.091,0.057,0.034,
                          0.019,0.01,0.005,0.003,
                          0.001,0.001)

  estimates <- smooth_deconvolve_estimate(incidence_data = toy_incidence_data,
                                          smoothing_method = "LOESS",
                                          deconvolution_method = "Richardson-Lucy delay distribution",
                                          estimation_method = "EpiEstim sliding window",
                                          delay_incubation = delay_distribution,
                                          estimation_window = 5,
                                          mean_serial_interval = 8,
                                          std_serial_interval  = 3,
                                          minimum_cumul_incidence = 0,
                                          block_size = 3,
                                          degree = 1,
                                          mean_Re_prior = 2.5,
                                          output_Re_only = FALSE,
                                          index_col = "date_index")

  reference_R_values <- c(NA,NA,NA,NA,NA,NA,NA,9.493,7.178,
                          5.892,5.114,4.591,4.201,3.887,3.609,
                          3.345,3.099,2.868,2.645,2.436,2.244,
                          2.063,1.895,1.74,1.593,1.454,1.324,
                          1.203,1.092,0.994,0.909,0.836,0.775,
                          0.722,0.674,0.629,0.586,NA,NA,NA,NA,NA)

  reference_indices <- -5:36

  expect_equal(estimates$R_mean, reference_R_values, tolerance = 1E-2)
  expect_equal(estimates$date_index, reference_indices)
})


#TODO skip on CRAN as it can fail by chance
test_that("get_block_bootstrapped_estimate yields consistent results on a toy example", {

  toy_incidence_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                          138,162,189,218,245,268,292,311,322,330,
                          332,324,312,297,276,256,236,214,192,170,
                          145,118,91,66)

  shape_incubation <-  2
  scale_incubation <- 1.2
  delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  estimates <- get_block_bootstrapped_estimate(incidence_data = toy_incidence_data,
                                               N_bootstrap_replicates = 100,
                                               smoothing_method = "LOESS",
                                               deconvolution_method = "Richardson-Lucy delay distribution",
                                               estimation_method = "EpiEstim sliding window",
                                               uncertainty_summary_method = "bagged mean - CI from bootstrap estimates",
                                               delay_incubation = delay_incubation,
                                               delay_onset_to_report = delay_onset_to_report,
                                               estimation_window = 3,
                                               mean_serial_interval = 4.8,
                                               std_serial_interval  = 2.3,
                                               minimum_cumul_incidence = 10,
                                               mean_Re_prior = 1,
                                               ref_date = as.Date("2020-02-04"),
                                               time_step = "day")

  reference_dates <- seq.Date(from = as.Date("2020-02-02"),
                              to = as.Date("2020-03-05"),
                              by="day")

  reference_R_mean_values <- c(4.981,3.595,3.027,2.766,2.619,2.508,2.408,
                               2.306,2.198,2.09,1.986,1.884,1.785,1.693,
                               1.604,1.518,1.437,1.362,1.29,1.22,1.153,
                               1.088,1.028,0.973,0.926,0.888,0.859,0.837,
                               0.819,0.804,0.787,0.769,0.75)

  reference_CI_down_values <- c(4.66,3.297,2.74,2.488,2.365,2.294,2.231,
                                2.149,2.051,1.95,1.854,1.767,1.687,1.61,
                                1.535,1.46,1.384,1.311,1.241,1.172,1.104,
                                1.039,0.975,0.912,0.854,0.803,0.762,0.727,
                                0.697,0.667,0.634,0.597,0.556)

  reference_CI_up_values <- c(5.053,3.666,3.096,2.832,2.676,2.552,2.449,
                              2.359,2.27,2.182,2.088,1.982,1.873,1.769,
                              1.669,1.575,1.49,1.413,1.339,1.265,1.189,
                              1.114,1.043,0.979,0.923,0.879,0.846,0.82,
                              0.799,0.779,0.758,0.733,0.708)


  expect_equal(estimates$date, reference_dates)
  expect_equal(estimates$Re_estimate, reference_R_mean_values, tolerance = 1E-1)
  expect_equal(estimates$CI_down, reference_CI_down_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up, reference_CI_up_values, tolerance = 1E-1)

  estimates <- get_block_bootstrapped_estimate(incidence_data = toy_incidence_data,
                                               N_bootstrap_replicates = 100,
                                               smoothing_method = "LOESS",
                                               deconvolution_method = "Richardson-Lucy delay distribution",
                                               estimation_method = "EpiEstim sliding window",
                                               uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                               delay_incubation = delay_incubation,
                                               delay_onset_to_report = delay_onset_to_report,
                                               estimation_window = 3,
                                               minimum_cumul_incidence = 0,
                                               mean_serial_interval = 4.8,
                                               std_serial_interval  = 2.3,
                                               mean_Re_prior = 1)

  reference_indices <- -2:30

  reference_R_original_values <- c(4.856,3.482,2.918,2.66,2.52,2.423,2.34,2.254,
                                   2.161,2.066,1.971,1.875,1.78,1.69,1.602,1.517,
                                   1.437,1.362,1.29,1.218,1.147,1.076,1.009,0.945,
                                   0.889,0.841,0.804,0.774,0.748,0.723,0.696,
                                   0.665,0.632)

  reference_CI_down_values <- c(4.682,3.321,2.769,2.523,2.401,2.327,2.264,2.19,
                                2.103,2.011,1.919,1.824,1.73,1.639,1.55,1.464,
                                1.384,1.311,1.243,1.177,1.112,1.049,0.989,0.932,
                                0.879,0.834,0.797,0.765,0.738,0.711,0.682,0.648,0.612)

  reference_CI_up_values <- c(5.031,3.642,3.067,2.797,2.639,2.519,2.415,2.318,2.218,
                              2.121,2.024,1.925,1.83,1.74,1.654,1.57,1.491,1.414,
                              1.338,1.26,1.181,1.103,1.028,0.959,0.898,0.849,0.811,
                              0.782,0.757,0.734,0.71,0.682,0.651)

  expect_equal(estimates$idx, reference_indices)
  expect_equal(estimates$Re_estimate, reference_R_original_values, tolerance = 1E-1)
  expect_equal(estimates$CI_down, reference_CI_down_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up, reference_CI_up_values, tolerance = 1E-1)
})

#TODO skip on CRAN as it can fail by chance
test_that("get_block_bootstrapped_estimate passes '...' arguments to inner functions properly", {

  toy_incidence_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                          138,162,189,218,245,268,292,311,322,330,
                          332,324,312,297,276,256,236,214,192,170,
                          145,118,91,66)

  shape_incubation <-  2
  scale_incubation <- 1.2
  delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  estimates <- get_block_bootstrapped_estimate(incidence_data = toy_incidence_data,
                                               N_bootstrap_replicates = 100,
                                               smoothing_method = "LOESS",
                                               deconvolution_method = "Richardson-Lucy delay distribution",
                                               estimation_method = "EpiEstim sliding window",
                                               uncertainty_summary_method = "bagged mean - CI from bootstrap estimates",
                                               delay_incubation = delay_incubation,
                                               delay_onset_to_report = delay_onset_to_report,
                                               estimation_window = 5,
                                               mean_serial_interval = 4.8,
                                               std_serial_interval  = 2.3,
                                               minimum_cumul_incidence = 0,
                                               block_size = 8,
                                               degree = 2,
                                               ref_date = as.Date("2020-02-04"),
                                               time_step = "day")

  reference_R_mean_values <- c(12.24,8.9,6.81,5.5,4.6,3.95,3.49,3.15,2.88,
                               2.66,2.48,2.31,2.15,2.02,1.89,1.76,1.64,
                               1.53,1.42,1.32,1.22,1.13,1.05,0.98,0.92,
                               0.86,0.8,0.74,0.67,0.6,0.51,0.42)

  reference_CI_down_values <- c(10.23,7.64,6.02,5.01,4.28,3.72,
                                3.28,2.95,2.69,2.49,2.34,2.2,2.06,
                                1.94,1.83,1.71,1.59,1.49,1.38,1.28,
                                1.19,1.1,1.03,0.96,0.9,0.84,0.78,
                                0.72,0.66,0.58,0.5,0.4)

  reference_CI_up_values <- c(14.26,10.16,7.6,5.99,4.92,
                              4.18,3.69,3.36,3.07,2.82,
                              2.62,2.43,2.25,2.09,1.95,
                              1.81,1.69,1.57,1.46,1.36,
                              1.26,1.17,1.08,1,0.93,0.87,
                              0.81,0.75,0.68,0.61,0.53,0.44)


  expect_equal(estimates$Re_estimate, reference_R_mean_values, tolerance = 1E-1)
  expect_equal(estimates$CI_down, reference_CI_down_values, tolerance = 1E-1)
  expect_equal(estimates$CI_up, reference_CI_up_values, tolerance = 1E-1)
})


test_that("get_infections_from_incidence handles partially-delayed data correctly",{
  toy_onset_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                      138,162,189,218,245,268,292,311,322,330,
                      332,324,312,297,276,256,236,214,192,170,
                      145,118,91,66,45,32,11,5,4,0,1,2,0)

  shape_incubation <-  2
  scale_incubation <- 1.2
  delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 3
  scale_onset_to_report <- 1.3
  delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  result_deconvolution <- get_infections_from_incidence(incidence_data = toy_onset_data,
                                smoothing_method = "LOESS",
                                deconvolution_method = "Richardson-Lucy delay distribution",
                                delay_incubation = delay_incubation,
                                is_partially_reported_data = TRUE,
                                delay_distribution_final_report = delay_onset_to_report,
                                output_infection_incidence_only = FALSE,
                                data_points_incl = 21,
                                degree = 1)

  reference_deconvolved_incidence <- c(16.417,19.763,22.179,26.041,32.397,40.832,
                               50.92,62.293,75.173,89.947,106.529,124.207,
                               142.552,161.947,182.9,203.07,221.026,239.061,
                               255.931,268.766,278.113,285.134,289.532,289.975,
                               286.109,278.965,269.38,256.736,240.254,221.707,
                               202.645,181.802,160.793,141.705,123.465,104.689,
                               85.982,68.065,49.84,31.857,15.892,3.778,
                               0,0,0,NA,NA,NA)

  expect_equal(result_deconvolution$deconvolved_incidence, reference_deconvolved_incidence, tolerance = 1E-1)
})

test_that("estimate_from_combined_observations returns consistent results",{
  toy_onset_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                      138,162,189,218,245,268,292,311,322,330,
                      332,324,312,297,276,256,236,214,192,170,
                      145,118,91,66,45,32,11,5,4,0,1,2,0)

  toy_case_confirmation_data <- c(11,12,21,23,2,14,49,61,65,45,66,45,
                                  40,8,61,38,1,3,45,66,12,52,27,3,54,
                                  10,18,54,12,48,67,62,54,3,29,10,52,
                                  61,33,39,55,8,64,51,65,34)

  shape_incubation <-  1
  scale_incubation <- 1.2
  delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 4
  scale_onset_to_report <- 2.3
  delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  results_estimation <- estimate_from_combined_observations(partially_delayed_incidence = toy_onset_data,
                                      fully_delayed_incidence = toy_case_confirmation_data,
                                      smoothing_method = "LOESS",
                                      deconvolution_method = "Richardson-Lucy delay distribution",
                                      estimation_method = "EpiEstim sliding window",
                                      delay_until_partial = delay_incubation,
                                      delay_from_partial_to_full = delay_onset_to_report,
                                      partial_observation_requires_full_observation = TRUE,
                                      ref_date = as.Date("2021-03-24"),
                                      time_step = "day",
                                      minimum_cumul_incidence = 0,
                                      output_Re_only = FALSE,
                                      data_points_incl = 21,
                                      degree = 1)



  reference_deconvolved_incidence <- c(59.2,61.7,64.2,69.8,76.9,84.4,92.1,101.5,112.9,126.2,
                                       140.8,156.9,175.5,195.1,214.5,235.3,255.6,272.3,287,
                                       301.6,314.4,323.1,328.8,332.7,333.8,330.9,324.9,317.1,
                                       305.9,290.8,273.5,255.4,237.4,219.6,201.4,183.3,165.1,
                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)

  reference_R_values <- c(NA,NA,NA,NA,3.5,2.4,2,1.7,1.6,1.6,1.6,1.6,1.6,1.6,
                          1.6,1.6,1.6,1.5,1.5,1.4,1.4,1.3,1.3,1.2,1.2,1.1,
                          1.1,1,1,0.9,0.9,0.9,0.8,0.8,0.7,0.7,0.7,
                          NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)

  expect_equal(results_estimation$combined_deconvolved_incidence, reference_deconvolved_incidence, tolerance = 1E-1)
  expect_equal(results_estimation$R_mean, reference_R_values, tolerance = 1E-1)
})

#TODO switch off on CRAN as can fail by chance
test_that("get_bootstrapped_estimate_from_combined_observations returns consistent results",{

  toy_onset_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                      138,162,189,218,245,268,292,311,322,330,
                      332,324,312,297,276,256,236,214,192,170,
                      145,118,91,66,45,32,11,5,4,0,1,2,0)

  toy_case_confirmation_data <- c(11,12,21,23,2,14,49,61,65,45,66,45,
                                  40,8,61,38,1,3,45,66,12,52,27,3,54,
                                  10,18,54,12,48,67,62,54,3,29,10,52,
                                  61,33,39,55,8,64,51,65,34)

  shape_incubation <-  1
  scale_incubation <- 1.2
  delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

  shape_onset_to_report <- 4
  scale_onset_to_report <- 2.3
  delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

  results_estimation <- get_bootstrapped_estimates_from_combined_observations(partially_delayed_incidence = toy_onset_data,
                                                            fully_delayed_incidence = toy_case_confirmation_data,
                                                            smoothing_method = "LOESS",
                                                            deconvolution_method = "Richardson-Lucy delay distribution",
                                                            estimation_method = "EpiEstim sliding window",
                                                            bootstrapping_method = "non-parametric block boostrap",
                                                            uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                                            N_bootstrap_replicates = 100,
                                                            delay_until_partial = delay_incubation,
                                                            delay_from_partial_to_full = delay_onset_to_report,
                                                            partial_observation_requires_full_observation = TRUE,
                                                            ref_date = as.Date("2021-03-24"),
                                                            time_step = "day",
                                                            output_Re_only = FALSE,
                                                            minimum_cumul_incidence = 0,
                                                            data_points_incl = 21,
                                                            degree = 1)


  reference_R_values <- c(3.5,2.4,2,1.7,1.6,1.6,1.6,1.6,1.6,1.6,
                          1.6,1.6,1.6,1.5,1.5,1.4,1.4,1.3,1.3,1.2,1.2,1.1,
                          1.1,1,1,0.9,0.9,0.9,0.8,0.8,0.7,0.7,0.7)

  reference_CI_up <- c(3.71,2.61,2.17,1.95,1.85,1.81,1.79,1.79,1.79,
                       1.78,1.76,1.73,1.68,1.62,1.56,1.49,1.43,1.38,
                       1.32,1.26,1.21,1.15,1.11,1.06,1.02,0.98,0.94,
                       0.9,0.86,0.83,0.81,0.79,0.77)

  reference_CI_down <- c(3.26,2.2,1.75,1.54,1.45,1.42,1.43,1.44,1.46,
                         1.48,1.49,1.49,1.48,1.46,1.41,1.36,1.31,1.25,
                         1.21,1.16,1.12,1.07,1.03,0.99,0.94,0.89,0.85,
                         0.8,0.76,0.72,0.69,0.65,0.61)

  expect_equal(results_estimation$Re_estimate, reference_R_values, tolerance = 1E-1)
  expect_equal(results_estimation$CI_up, reference_CI_up, tolerance = 1E-1)
  expect_equal(results_estimation$CI_down, reference_CI_down, tolerance = 1E-1)
})

test_that("get_bootstrapped_estimate_from_combined_observations can deal with empirical delay data",{

  ref_date <- as.Date("2021-03-24")
  n_days <- 50
  shape_initial_delay <- 6
  scale_initial_delay <- 1.5
  distribution_initial_delay <- list(name = "gamma", shape = shape_initial_delay, scale = scale_initial_delay)
  seed <- 7543265

  generated_empirical_delays <- .generate_delay_data(origin_date = ref_date,
                                                     n_time_steps = n_days,
                                                     ratio_delay_end_to_start = 1.5,
                                                     distribution_initial_delay = distribution_initial_delay,
                                                     seed = seed)

  toy_onset_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                      138,162,189,218,245,268,292,311,322,330,
                      332,324,312,297,276,256,236,214,192,170,
                      145,118,91,66,45,32,11,5,4,0,1,2,0)

  toy_case_confirmation_data <- c(11,12,21,23,2,14,49,61,65,45,66,45,
                                  40,8,61,38,1,3,45,66,12,52,27,3,54,
                                  10,18,54,12,48,67,62,54,3,29,10,52,
                                  61,33,39,55,8,64,51,65,34)

  shape_incubation <-  2
  scale_incubation <- 1.2
  delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

  results_estimation <- get_bootstrapped_estimates_from_combined_observations(partially_delayed_incidence = toy_onset_data,
                                                                             fully_delayed_incidence = toy_case_confirmation_data,
                                                                             smoothing_method = "LOESS",
                                                                             deconvolution_method = "Richardson-Lucy delay distribution",
                                                                             estimation_method = "EpiEstim sliding window",
                                                                             bootstrapping_method = "non-parametric block boostrap",
                                                                             uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                                                             N_bootstrap_replicates = 20, # to speed things up
                                                                             delay_until_partial = delay_incubation,
                                                                             delay_from_partial_to_full = generated_empirical_delays,
                                                                             partial_observation_requires_full_observation = TRUE,
                                                                             ref_date = ref_date,
                                                                             time_step = "day",
                                                                             output_Re_only = FALSE,
                                                                             minimum_cumul_incidence = 0,
                                                                             data_points_incl = 21,
                                                                             degree = 1)


  reference_R_values <- c(3.54,2.43,2.01,1.86,1.81,1.79,1.78,
                          1.77,1.76,1.75,1.71,1.65,1.58,1.52,
                          1.45,1.39,1.33,1.28,1.23,1.17,1.12,
                          1.07,1.02,0.97,0.93,0.88,0.84,0.8,
                          0.76,0.72,0.69)

  expect_equal(results_estimation$Re_estimate, reference_R_values, tolerance = 1E-1)
})


