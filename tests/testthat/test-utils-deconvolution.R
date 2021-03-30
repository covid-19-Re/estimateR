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





