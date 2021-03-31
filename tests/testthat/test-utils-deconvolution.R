#TODO TEST that:
# 1) build_delay_distribution throws error when unsupported distribution_type is thrown in
#and when unsuitable parameter values are thrown in (not numeric, or negative values for instance)
# 2) get_matrix_empirical_waiting_time_distr  reports consistent results on a toy example
#

test_that("build_delay_distribution returns a vector whose elements sum up to 1", {
  N <- 100
  shapes <- stats::runif(N, min = 0, max = 10)
  scales <- stats::runif(N, min = 0, max = 10)

  delay_distribution_vectors <- sapply(1:N, function(x) {
                                  build_delay_distribution(parm1 = shapes[x],
                                                          parm2 = scales[x],
                                                          distribution_type = "gamma",
                                                          max_quantile = 0.9999)
    })

  max_difference_to_1 <- max(abs(sapply(delay_distribution_vectors, sum) - 1))

  expect_equal(max_difference_to_1, 0, tolerance = 1E-4)
})


test_that(".convolve_delay_distribution_vectors returns a vector whose elements sum up to 1", {
  N <- 100
  shapes <- stats::runif(2*N, min = 0, max = 10)
  scales <- stats::runif(2*N, min = 0, max = 10)

  delay_distribution_vectors <- sapply(1:(2*N), function(x) {
    build_delay_distribution(parm1 = shapes[x],
                             parm2 = scales[x],
                             distribution_type = "gamma",
                             max_quantile = 0.9999)
  })

  convolved_distribution_vectors <- sapply(1:N, function(x) {
    .convolve_delay_distribution_vectors(delay_distribution_vectors[[2*x]],
                                         delay_distribution_vectors[[2*x-1]])
  })

  max_difference_to_1 <- max(abs(sapply(delay_distribution_vectors, sum) - 1))

  expect_equal(max_difference_to_1, 0, tolerance = 1E-4)
})


test_that(".convolve_delay_distribution_vectors returns correct output on a simple example", {

  delay_distribution_vector_1 <- c(0, 0.25, 0.1, 0.65)
  delay_distribution_vector_2 <- c(0.2, 0.2, 0.3, 0.3)
  ref_convolved_output <- c(0, 0.05, 0.07, 0.225, 0.235, 0.225, 0.195, 0)

  convolved_output <- .convolve_delay_distribution_vectors(delay_distribution_vector_1,
                                       delay_distribution_vector_2)

  expect_equal(convolved_output, ref_convolved_output, tolerance = 1E-4)
})

test_that(".combine_incubation_with_reporting_delay returns same output as empirical method of convoluting gammas", {

  shape_incubation = 3.2
  scale_incubation = 1.3

  shape_onset_to_report = 2.7
  scale_onset_to_report = 1.6

  convolved_output <- combine_incubation_with_reporting_delay(parm1_incubation=shape_incubation,
                                                              parm2_incubation=scale_incubation,
                                                              parm1_onset_to_report=shape_onset_to_report,
                                                              parm2_onset_to_report=scale_onset_to_report,
                                                              distribution_type_incubation="gamma",
                                                              distribution_type_onset_to_report="gamma",
                                                              max_quantile = 0.9999)

  empirical_convolution_result <- c(0,9e-04,0.00947,0.03214,0.06438,0.09523,0.11566,
                                    0.12255,0.11742,0.1043,0.08721,0.06951,0.05329,
                                    0.03951,0.02841,0.01991,0.01371,0.00925,0.00616,
                                    0.00402,0.0026,0.00165,0.00105,0.00065,0.00039,
                                    0.00025,0.00015,9e-05,6e-05,3e-05,2e-05,1e-05,
                                    1e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

 padded_convolved_output <- c(convolved_output, rep(0, times = length(empirical_convolution_result) - length(convolved_output)))

 absolute_diff <- abs(padded_convolved_output - empirical_convolution_result)

 expect_equal(max(absolute_diff), 0, tolerance = 1E-3)
})





