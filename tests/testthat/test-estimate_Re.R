test_that("estimate_Re yields consistent results on a toy example", {
  incidence_data <- c(
    1, 2, 12, 32, 34, 45, 87, 134, 230, 234, 222, 210, 190, 259,
    351, 453, 593, 603, 407, 348, 304, 292, 256, 229,
    132, 98, 86, 54, 39, 23, 3, 2, 12, 14
  )

  estimated_Re <- estimate_Re(
    incidence_data = incidence_data,
    estimation_method = "EpiEstim sliding window",
    minimum_cumul_incidence = 0,
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    mean_Re_prior = 1
  )

  reference_values <- c(
    23.12, 10.59, 7.01, 6.18, 6.4, 5.32, 3.83, 2.46, 1.67, 1.42, 1.51, 1.84,
    2.22, 2.31, 1.89, 1.32, 0.89, 0.73, 0.65, 0.62, 0.53, 0.43, 0.33, 0.29,
    0.26, 0.2, 0.14, 0.08, 0.06, 0.14
  )
  reference_offset <- 4

  expect_equal(.get_values(estimated_Re), reference_values, tolerance = 1E-2)
  expect_identical(.get_offset(estimated_Re), reference_offset)
})

test_that("estimate_Re is consistent with piecewise estimates", {
  incidence_data <- c(
    1, 2, 12, 32, 34, 45, 87, 134, 230, 234, 222, 210, 190, 259,
    351, 453, 593, 603, 407, 348, 304, 292, 256, 229,
    132, 98, 86, 54, 39, 23, 3, 2, 12, 14
  )

  piecewise_estimated_Re <- estimate_Re(
    incidence_data = incidence_data,
    estimation_method = "EpiEstim piecewise constant",
    minimum_cumul_incidence = 0,
    interval_length = 7,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    mean_Re_prior = 1
  )

  reference_piecewise_values <- c(
    7.46, 7.46, 7.46, 7.46, 7.46, 7.46, 7.46,
    2.03, 2.03, 2.03, 2.03, 2.03, 2.03, 2.03,
    1.28, 1.28, 1.28, 1.28, 1.28, 1.28, 1.28,
    0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41,
    0.12, 0.12, 0.12, 0.12, 0.12
  )

  reference_piecewise_offset <- 1

  expect_equal(.get_values(piecewise_estimated_Re), reference_piecewise_values, tolerance = 1E-2)
  expect_identical(.get_offset(piecewise_estimated_Re), reference_piecewise_offset)

  piecewise_estimated_Re <- estimate_Re(
    incidence_data = incidence_data,
    estimation_method = "EpiEstim piecewise constant",
    minimum_cumul_incidence = 0,
    interval_ends = c(9, -2, 0, 1, 14, 49, 78, 34),
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    mean_Re_prior = 1
  )

  reference_piecewise_values <- c(
    7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1,
    1.83, 1.83, 1.83, 1.83, 1.83, 0.83, 0.83,
    0.83, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83,
    0.83, 0.83, 0.83, 0.83, 0.83, 0.83, 0.83,
    0.83, 0.83, 0.83, 0.83
  )

  reference_piecewise_offset <- 1

  expect_equal(.get_values(piecewise_estimated_Re), reference_piecewise_values, tolerance = 1E-2)
  expect_identical(.get_offset(piecewise_estimated_Re), reference_piecewise_offset)
})

test_that("estimate_Re outputs HPD intervals when asked", {
  incidence_data <- c(
    1, 2, 12, 32, 34, 45, 87, 134, 230, 234, 222, 210, 190, 259,
    351, 453, 593, 603, 407, 348, 304, 292, 256, 229,
    132, 98, 86, 54, 39, 23, 3, 2, 12, 14
  )

  estimated_Re <- estimate_Re(
    incidence_data = incidence_data,
    estimation_method = "EpiEstim sliding window",
    minimum_cumul_incidence = 0,
    estimation_window = 3,
    mean_serial_interval = 4.8,
    std_serial_interval = 2.3,
    mean_Re_prior = 1,
    output_HPD = TRUE,
    simplify_output = TRUE
  )

  reference_index <- c(
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22,
    23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33
  )

  reference_values_lowHPD <- c(
    18.28, 8.71, 5.99, 5.46, 5.82,
    4.91, 3.55, 2.28, 1.54, 1.31, 1.41,
    1.73, 2.11, 2.2, 1.79, 1.25, 0.83,
    0.68, 0.61, 0.58, 0.49, 0.39, 0.3,
    0.25, 0.22, 0.17, 0.11, 0.05, 0.04, 0.09
  )

  expect_identical(estimated_Re$idx, reference_index)
  expect_equal(estimated_Re$Re_lowHPD, reference_values_lowHPD, tolerance = 1E-2)
})
