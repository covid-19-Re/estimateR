test_that("smooth_deconvolve_estimate yields consistent results on a toy example", {

  toy_incidence_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                          138,162,189,218,245,268,292,311,322,330,
                          332,324,312,297,276,256,236,214,192,170,
                          145,118,91,66)

  delay_distribution <- c(0,0.015,0.09,0.168,
                          0.195,0.176,0.135,0.091,0.057,0.034,
                          0.019,0.01,0.005,0.003,
                          0.001,0.001)

  estimates <- smooth_deconvolve_estimate(toy_incidence_data,
                                         delay_distribution,
                                         smoothing_method = "LOESS",
                                         deconvolution_method = "Richardson-Lucy delay distribution",
                                         estimation_method = "EpiEstim sliding window",
                                         estimation_window = 3,
                                         mean_serial_interval = 4.8,
                                         std_serial_interval  = 2.3,
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



test_that("get_block_bootstrapped_estimate yields consistent results on a toy example", {

  toy_incidence_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                          138,162,189,218,245,268,292,311,322,330,
                          332,324,312,297,276,256,236,214,192,170,
                          145,118,91,66)

  estimates <- get_block_bootstrapped_estimate(toy_incidence_data,
                                               N_bootstrap_replicates = 100,
                                          smoothing_method = "LOESS",
                                          deconvolution_method = "Richardson-Lucy delay distribution",
                                          estimation_method = "EpiEstim sliding window",
                                          shape_incubation = 2,
                                          scale_incubation = 1.2,
                                          shape_onset_to_report = 3,
                                          scale_onset_to_report = 1.3,
                                          estimation_window = 3,
                                          mean_serial_interval = 4.8,
                                          std_serial_interval  = 2.3,
                                          mean_Re_prior = 1,
                                          ref_date = as.Date("2020-02-04"),
                                          time_step = "day")

  estimated_R_mean <- estimates %>%
    dplyr::select(date, R_mean) %>%
    na.omit %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(mean_R = mean(R_mean)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(mean_R)

  reference_R_values <- c(4.981,3.595,3.027,2.766,2.619,2.508,2.408,
                          2.306,2.198,2.09,1.986,1.884,1.785,1.693,
                          1.604,1.518,1.437,1.362,1.29,1.22,1.153,
                          1.088,1.028,0.973,0.926,0.888,0.859,0.837,
                          0.819,0.804,0.787,0.769,0.75)
  expect_equal(estimated_R_mean, reference_R_values, tolerance = 1E-1)
})
