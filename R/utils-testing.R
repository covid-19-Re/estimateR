#' Utility function to print vectors in a copy-pastable format
#'
#' @param a numeric vector
#' @param digits integer. Number of digits to round.
#'
#' @return Print vector as string
.print_vector <- function(a, digits = 2) {
  cat(paste0("c(", paste(round(a, digits = digits), collapse = ","), ")"))
}

#' Generate artificial delay data.
#'
#' This utility can be used to build toy examples to test functions dealing with empirical delay data.
#' It is very basic in what it simulates.
#' A random walk is simulated over \code{n_time_steps}, representing the incidence through time.
#' The result of this simulation is offset so that all values are positive.
#' Then, for each time step, \code{n} samples from a delay distribution are taken,
#' with \code{n} being the incidence value at this time step.
#' The random draws are then multiplied by a factor (>1 or <1) to simulate
#' a gradual shift in the delay distribution through time.
#' This multiplication factor is calculated
#' by linearly interpolating between 1 (at the first time step),
#' and \code{delay_ratio_start_to_end} linearly,
#' from 1 at the first time step to \code{ratio_delay_end_to_start}
#' at the last time step.
#'
#' @param origin_date Date of first infection.
#' @param n_time_steps interger. Number of time steps to generate delays over
#' @param ratio_delay_end_to_start numeric value.
#' Shift in delay distribution from start to end.
#' @param distribution_initial_delay Distribution in list format.
#' @param seed integer. Optional RNG seed.
#' @inherit dating
#'
#' @return dataframe. Simulated delay data.
.generate_delay_data <- function(origin_date = as.Date("2020-02-01"),
                                 n_time_steps = 100,
                                 time_step = "day",
                                 ratio_delay_end_to_start = 2,
                                 distribution_initial_delay = list(name = "gamma", shape = 6, scale = 5),
                                 seed = NULL) {
  .are_valid_argument_values(list(
    list(origin_date, "date"),
    list(n_time_steps, "positive_integer"),
    list(time_step, "time_step"),
    list(ratio_delay_end_to_start, "number"),
    list(distribution_initial_delay, "distribution"),
    list(seed, "null_or_int")
  ))

  if (!is.null(seed)) {
    set.seed(seed)
  }

  random_pop_size_walk <- cumsum(sample(c(-1, 0, 1), n_time_steps, TRUE))
  pop_size <- random_pop_size_walk - min(0, min(random_pop_size_walk)) + 1

  event_dates <- seq.Date(from = origin_date, length.out = n_time_steps, by = time_step)

  delays <- lapply(1:n_time_steps, function(i) {
    # Sample a number of draws from the gamma delay distribution, based on the pop size at time step i
    raw_sampled_delays <- .sample_from_distribution(
      distribution = distribution_initial_delay,
      n = pop_size[i]
    )
    # Multiply these samples by a factor that accounts for the linear inflation or deflation of delays.
    sampled_delays <- round(raw_sampled_delays * (1 + i / n_time_steps * (ratio_delay_end_to_start - 1)))

    return(tibble::tibble(event_date = event_dates[i], report_delay = sampled_delays))
  })

  return(dplyr::bind_rows(delays))
}

#' Utility function that generates delay data, assuming a different delay between event and observation for each individual day.
#' It then generates a delay matrix and computes the RMSE between the parameters of the gamma distributions passed as arguments and the ones recovered from the delay matrix.
#' The shapes and scales of the gamma distributions are specified as parameters, and the number of timesteps is assumed to be equal to the length of these vectors.
#'
#' This funciotn is useful for testing purposes.
#' @param original_distribution_shapes vector. Specifies the shapes for the gamma distributions.
#' @param original_distribution_scales vector. Specifies the scales for the gamma distributions.
#' @param nr_distribution_samples integer. How many cases to be sampled for each timestep.
#'
#' @return A list with the computed RMSE. It has two elements: $shape_rmse and $scale_rmse
.delay_distribution_matrix_rmse_compute <- function(original_distribution_shapes, original_distribution_scales, nr_distribution_samples = 500){

  #Create a vector with all dates in observation interval
  start_date <- as.Date('2021/04/01')
  time_steps <- length(original_distribution_shapes)
  end_date <- start_date + time_steps
  available_dates <- seq(start_date, end_date, by="day")

  #Build the delay data; Events on each individual day are assumed to be observed according to a different gamma distribution, as specified by original_distribution_shapes and original_distribution_scales,
  sampled_report_delays <- c()
  report_dates <- as.Date(c())
  for (i in 1:time_steps){
    new_sampled_report_delays <- .sample_from_distribution(list(name="gamma", shape=original_distribution_shapes[i], scale=original_distribution_scales[i]), nr_distribution_samples)
    sampled_report_delays <- c(sampled_report_delays, new_sampled_report_delays)
    new_report_dates <- rep(available_dates[i], nr_distribution_samples)
    report_dates <- c(report_dates, new_report_dates)
  }
  delay_data <- dplyr::tibble(event_date = report_dates, report_delay = sampled_report_delays)
  result <- get_matrix_from_empirical_delay_distr(delay_data, time_steps, fit = "gamma", return_fitted_distribution = TRUE)

  delay_matrix <- result$matrix
  distrib_list <- result$distributions

  #Get the shapes and scales of the gamma distributions fitted by the get_matrix_from_empirical_delay_distr function
  distribution_shapes <- c()
  distribution_scales <- c()

  for (distribution in distrib_list){
    distribution_shapes <- c(distribution_shapes, distribution$shape)
    distribution_scales <- c(distribution_scales, distribution$scale)
  }

  #Compute the RMSE between the desired gamma distribution shapes and scales, and the ones obtained by the get_matrix_from_empirical_delay_distr function
  start_index <- length(distribution_shapes) - length(original_distribution_shapes) + 1
  shape_rmse <- Metrics::rmse(distribution_shapes[start_index:length(distribution_shapes)], original_distribution_shapes)/mean(original_distribution_shapes)
  scale_rmse <- Metrics::rmse(distribution_scales[start_index:length(distribution_scales)], original_distribution_scales)/mean(original_distribution_scales)

  return(list(shape_rmse=shape_rmse, scale_rmse=scale_rmse))
}

#' Test whether a delay matrix columns sum to less than one.
#'
#' @param matrix input matrix
#' @param full_cols number of full columns.
#' Full columns must sum to 1.
#' @param tolerance tolerance on the result
#'
#' @return TRUE if test passed
expect_delay_matrix_sums_lte_1 <- function(matrix, full_cols = 0, tolerance = 1E-3) {
  if (full_cols > 0 && full_cols < ncol(matrix)) {
    sums_full_cols <- apply(matrix[, 1:full_cols], MARGIN = 2, FUN = sum)
    testthat::expect_equal(sums_full_cols, rep(1, times = length(sums_full_cols)), tolerance = tolerance)
  }

  sum_all_cols <- apply(matrix, MARGIN = 2, FUN = sum)
  testthat::expect_lte(max(abs(sum_all_cols)), 1)
}
