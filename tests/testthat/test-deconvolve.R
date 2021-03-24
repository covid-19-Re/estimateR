test_that("deconvolve_incidence yields consistent results on a toy example", {
  toy_incidence_data <- c(6,8,10,13,17,22,31,41,52,65,80,97,116,
                          138,162,189,218,245,268,292,311,322,330,
                          332,324,312,297,276,256,236,214,192,170,
                          145,118,91,66)

  delay_distribution <- c(0,0.015,0.09,0.168,
                          0.195,0.176,0.135,0.091,0.057,0.034,
                          0.019,0.01,0.005,0.003,
                          0.001,0.001)

  deconvolved_incidence <- deconvolve_incidence(incidence_data = toy_incidence_data,
                                                deconvolution_method = "Richardson-Lucy delay distribution",
                                                delay_distribution = delay_distribution,
                                                time_units_in_the_past = 30,
                                                threshold_chi_squared = 1,
                                                max_iterations = 100)

  reference_values <- c(6,8,9,12,15,19,27,37,47,60,75,92,111,133,
                        158,186,217,247,272,299,321,334,343,345,337,
                        323,306,283,261,239,216,193,169,143,114,85,59)
  reference_offset <- -5

  expect_equal(.get_values(deconvolved_incidence), reference_values, tolerance = 1)
  expect_identical(.get_offset(deconvolved_incidence), reference_offset)
})
