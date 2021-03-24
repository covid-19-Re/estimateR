test_that("estimate_Re yields consistent results on a toy example", {
  incidence_data <- c(1,2,12,32,34,45,87,134, 230,234,222,210, 190, 259,
                      351, 453, 593, 603,407, 348, 304, 292, 256, 229,
                      132, 98, 86, 54, 39, 23, 3,2,12,14)

  estimated_Re <- estimate_Re(incidence_data = incidence_data,
              estimation_method = "EpiEstim sliding window",
              minimul_cumul_incidence = 0,
              estimation_window = 3,
              mean_serial_interval = 4.8,
              std_serial_interval  = 2.3,
              mean_Re_prior = 1)

  reference_values <- c(23.12,10.59,7.01,6.18,6.4,5.32,3.83,2.46,1.67,1.42,1.51,1.84,
                        2.22,2.31,1.89,1.32,0.89,0.73,0.65,0.62,0.53,0.43,0.33,0.29,
                        0.26,0.2,0.14,0.08,0.06,0.14)
  reference_offset <- 4

  expect_equal(.get_values(estimated_Re), reference_values, tolerance = 1E-2)
  expect_identical(.get_offset(estimated_Re), reference_offset)
})
