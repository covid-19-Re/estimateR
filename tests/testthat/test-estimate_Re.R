#TODO list what needs to be tested
#TEST that keeps reference value

test_that("multiplication works", {
  expect_equal(2 * 2, 4)

  #TODO continue

  incidence_data <- c(1,2,12,32,34,45,87,134, 230,234,222,210, 190, 259,
                      351, 453, 593, 603,407, 348, 304, 292, 256, 229,
                      132, 98, 86, 54, 39, 23, 3,2,12,14)

  input <- .get_module_input(incidence_data)

  .estimate_Re_EpiEstim_sliding_window(input,
                                       minimul_cumul_incidence = 0,
                                       estimation_window = 3,
                                       mean_serial_interval = 4.8,
                                       std_serial_interval  = 2.3,
                                       mean_Re_prior = 1)


  a <- estimate_Re(incidence_data = incidence_data,
              estimation_method = "EpiEstim sliding window",
              minimul_cumul_incidence = 0,
              estimation_window = 3,
              mean_serial_interval = 4.8,
              std_serial_interval  = 2.3,
              mean_Re_prior = 1)
})
