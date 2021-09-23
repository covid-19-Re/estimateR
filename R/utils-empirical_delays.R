#' Get delay entries corresponding to a round number of weeks (months...)
#'
#' The main use of this function is to allow comparison with legacy code.
#'
#' @inherit empirical_delay_data_format
#' @inherit delay_empirical
#'
#' @return A vector of length at least \code{min_number_cases}, containing
#' records of delays.
.get_delays_over_full_time_units <- function(delays,
                                             date_of_interest,
                                             num_steps_in_a_unit = 7,
                                             min_number_cases) {
  recent_counts <- delays %>%
    dplyr::arrange(dplyr::desc(.data$event_date)) %>%
    dplyr::filter(.data$event_date <= date_of_interest)

  if (nrow(recent_counts) < min_number_cases) {
    first_observation_dates <- delays %>%
      dplyr::arrange(.data$event_date) %>%
      dplyr::slice_head(n = min_number_cases) %>%
      dplyr::pull(.data$event_date)

    max_date <- max(first_observation_dates)
    num_steps_since_start <- trunc(as.double(max_date - min(delays$event_date), units = "auto"))

    if ((num_steps_since_start %% num_steps_in_a_unit) == 0) {
      max_date_with_full_weeks <- max_date
    } else {
      max_date_with_full_weeks <- max_date + num_steps_in_a_unit - (num_steps_since_start %% num_steps_in_a_unit)
    }

    recent_counts_distribution <- delays %>%
      dplyr::filter(.data$event_date <= max_date_with_full_weeks) %>%
      dplyr::pull(.data$report_delay)
  } else {
    first_observation_dates <- recent_counts %>%
      dplyr::slice_head(n = min_number_cases) %>%
      dplyr::pull(.data$event_date)

    min_date <- min(first_observation_dates)
    num_steps_since_start <- trunc(as.double(date_of_interest - min_date, units = "auto"))

    if ((num_steps_since_start %% num_steps_in_a_unit) == 0) {
      min_date_with_full_weeks <- min_date
    } else {
      min_date_with_full_weeks <- min_date - num_steps_in_a_unit + (num_steps_since_start %% num_steps_in_a_unit)
    }
    recent_counts_distribution <- delays %>%
      dplyr::filter(
        .data$event_date >= min_date_with_full_weeks,
        .data$event_date <= date_of_interest
      ) %>%
      dplyr::pull(.data$report_delay)
  }
  return(recent_counts_distribution)
}


#' Build matrix of delay distributions through time from empirical delay data.
#'
#' This function takes a record of delays between events and their observations
#' and builds a discretized delay distribution matrix from this record.
#' The discretized delay distribution matrix
#' is required for the application of the Richardson-Lucy algorithm.
#' The main benefit of providing empirical delay data to an \code{estimateR} analysis,
#' as opposed to specifiying a delay as a single distribution
#' (whether a fitted or empirical distribution) is that the variability of
#' the delays through time is used to inform the analysis and provide more accurate estimates.
#' If the average of delays has shifted from 5 days to 3 days between the beginning and end
#' of epidemic of interest, this will be reflected in the recorded empirical delays
#' and will be accounted for by \code{estimateR} when estimating the reproductive number.
#'
#' The \code{ref_date} argument here will be understood as the date of the first record in the incidence data
#' for which the empirical delay data will be used.
#' If \code{ref_date} is not provided, the reference date will be taken as being the earliest date in the
#' \code{event_date} column of \code{empirical_delays}. In other words, the date of the first record in the incidence data
#' will be assumed to be the same as the date of the first record in the empirical delay data.
#' If this is not the case in your analysis, make sure to specify a \code{ref_date} argument.
#'
#' @example man/examples/get_matrix_from_empirical_delay_distr.R
#'
#' @inherit empirical_delay_data_format
#'
#' @inherit delay_empirical
#' @inheritParams dating
#' @inheritParams .get_delay_matrix_column
#'
#' @return a discretized delay distribution matrix.
#' @export
get_matrix_from_empirical_delay_distr <- function(empirical_delays,
                                                  n_report_time_steps,
                                                  ref_date = NULL,
                                                  time_step = "day",
                                                  min_number_cases = NULL,
                                                  upper_quantile_threshold = 0.99,
                                                  min_number_cases_fraction = 0.2,
                                                  min_min_number_cases = 500,
                                                  fit = "none",
                                                  return_fitted_distribution = FALSE,
                                                  num_steps_in_a_unit = NULL) {
  .are_valid_argument_values(list(
    list(empirical_delays, "empirical_delay_data"),
    list(n_report_time_steps, "positive_integer"),
    list(ref_date, "null_or_date"),
    list(time_step, "time_step"),
    list(min_number_cases, "null_or_int"),
    list(upper_quantile_threshold, "numeric_between_zero_one"),
    list(min_number_cases_fraction, "numeric_between_zero_one"),
    list(min_min_number_cases, "positive_integer"),
    list(fit, "delay_matrix_column_fit"),
    list(return_fitted_distribution, "boolean"),
    list(num_steps_in_a_unit, "null_or_int")
  ))


  if (is.null(ref_date)) {
    ref_date <- min(dplyr::pull(empirical_delays, .data$event_date), na.rm = TRUE)
  }

  all_report_dates <- seq.Date(from = ref_date, by = time_step, length.out = n_report_time_steps)

  # Ignore the delay data that is posterior to the last incidence report date.
  empirical_delays <- empirical_delays %>%
    dplyr::filter(.data$event_date <= max(all_report_dates))

  # Set the 'min_number_cases' parameter if not set by the user
  # TODO make this 'min_number_cases' depend on the length of the time_series.
  if (is.null(min_number_cases)) {
    min_number_cases <- min_number_cases_fraction * nrow(empirical_delays)
    min_number_cases <- max(min_number_cases, min_min_number_cases)
  }

  # Find the threshold for right-truncation
  # No time-variation beyond this threshold due to the fraction of unsampled individuals when nearing the last sampling date
  # TODO put the search for threshold_right_truncation in separate utility function
  delay_counts <- empirical_delays %>%
    dplyr::select(.data$report_delay) %>%
    dplyr::group_by(.data$report_delay) %>%
    dplyr::summarise(counts = dplyr::n(), .groups = "drop")

  threshold_right_truncation <- delay_counts %>%
    dplyr::mutate(cumul_freq = cumsum(.data$counts) / sum(.data$counts)) %>%
    dplyr::filter(.data$cumul_freq > upper_quantile_threshold) %>%
    utils::head(n = 1) %>%
    dplyr::pull(.data$report_delay)

  # We left-pad the time range with a number of time steps corresponding in the initial shift in the deconvolution.
  # TODO it may be simpler to just do the augmentation during the deconvolution step
  initial_shift <- ceiling(stats::quantile(empirical_delays$report_delay, probs = 0.99, na.rm = T))[1]

  # Left-pad the dates we are looking at to account for shift between event dates and observation dates.
  all_dates <- c(
    rev(seq.Date(from = ref_date, by = paste0("-1 ", time_step), length.out = initial_shift + 1)),
    seq.Date(from = ref_date, by = time_step, length.out = n_report_time_steps)[-1]
  )

  n_time_steps <- n_report_time_steps + initial_shift

  delay_distribution_matrix <- matrix(0, nrow = n_time_steps, ncol = n_time_steps)

  last_varying_col <- dplyr::if_else(n_time_steps > threshold_right_truncation, n_time_steps - threshold_right_truncation, n_time_steps)

  distrib_list <- list() # needed for the test that checks if get_matrix_from_empirical_delay_distr returns a matrix with the expected distributions when using fit = "gamma"

  # Shuffle rows so as to get rid of potential biases associated
  shuffled_delays <- empirical_delays %>%
    dplyr::slice(sample(1:dplyr::n()))

  # Populate the delay_distribution_matrix by column
  if (n_time_steps > threshold_right_truncation) {
    for (i in 1:last_varying_col) {
      if (is.null(num_steps_in_a_unit)) {
        # TODO take out in internal function to reduce duplication
        recent_counts <- shuffled_delays %>%
          dplyr::arrange(dplyr::desc(.data$event_date)) %>%
          dplyr::filter(.data$event_date <= all_dates[i])

        if (nrow(recent_counts) >= min_number_cases) {
          # If enough data points before date of interest,
          # take most recent observations before this date.

          recent_counts_distribution <- recent_counts %>%
            dplyr::slice_head(n = min_number_cases) %>%
            dplyr::pull(.data$report_delay)
        } else {
          # Otherwise, take 'min_number_of_cases' observations,
          # even after date of interest.
          recent_counts_distribution <- shuffled_delays %>%
            dplyr::arrange(.data$event_date) %>%
            dplyr::slice_head(n = min_number_cases) %>%
            dplyr::pull(.data$report_delay)
        }
      } else {
        recent_counts_distribution <- .get_delays_over_full_time_units(
          delays = shuffled_delays,
          date_of_interest = all_dates[i],
          num_steps_in_a_unit = num_steps_in_a_unit,
          min_number_cases = min_number_cases
        )
      }

      result <- .get_delay_matrix_column(recent_counts_distribution, fit, col_number = i, n_time_steps, return_fitted_distribution)
      if (is.list(result)) {
        distrib_list[[i]] <- result$distribution
        new_column <- result$column
      } else {
        new_column <- result
      }
      delay_distribution_matrix[, i] <- new_column
    }
  } else { # if n_time_steps <= threshold_right_truncation

    if (is.null(num_steps_in_a_unit)) {
      recent_counts <- shuffled_delays %>%
        dplyr::arrange(dplyr::desc(.data$event_date)) %>%
        dplyr::filter(.data$event_date <= all_dates[1])

      if (nrow(recent_counts) >= min_number_cases) {
        # If enough data points before date of interest,
        # take most recent observations before this date.

        recent_counts_distribution <- recent_counts %>%
          dplyr::slice_head(n = min_number_cases) %>%
          dplyr::pull(.data$report_delay)
      } else {
        # Otherwise, take 'min_number_of_cases' observations,
        # even after date of interest.
        recent_counts_distribution <- shuffled_delays %>%
          dplyr::arrange(.data$event_date) %>%
          dplyr::slice_head(n = min_number_cases) %>%
          dplyr::pull(.data$report_delay)
      }
    } else {
      recent_counts_distribution <- .get_delays_over_full_time_units(
        delays = shuffled_delays,
        date_of_interest = all_dates[1],
        num_steps_in_a_unit = num_steps_in_a_unit,
        min_number_cases = min_number_cases
      )
    }



    result <- .get_delay_matrix_column(recent_counts_distribution, fit, col_number = 1, n_time_steps, return_fitted_distribution)
    if (is.list(result)) {
      distrib_list[[i]] <- result$distribution
      new_column <- result$column
    } else {
      new_column <- result
    }

    for (i in 0:(last_varying_col - 1)) {
      delay_distribution_matrix[, i + 1] <- c(rep(0, times = i), new_column[1:(length(new_column) - i)])
      if (fit == "gamma") {
        distrib_list <- append(distrib_list, distrib_list[length(distrib_list)])
      }
    }
  }

  if (last_varying_col < n_time_steps) {
    for (j in 1:threshold_right_truncation) {
      delay_distribution_matrix[, i + j] <- c(rep(0, times = j), delay_distribution_matrix[1:(nrow(delay_distribution_matrix) - j), i])
      if (fit == "gamma") {
        distrib_list <- append(distrib_list, distrib_list[length(distrib_list)])
      }
    }
  }

  if (return_fitted_distribution) {
    return(list(matrix = delay_distribution_matrix, distributions = distrib_list))
  }

  return(delay_distribution_matrix)
}


#' Build a specific column of the delay distribution matrix
#'
#' @param recent_counts_distribution numeric vector of report delays, as used in \code{get_matrix_from_empirical_delay_distr}
#' @param fit string. Can be either "none" or "gamma". Specifies the type of fitting applied to the computed column
#' @param col_number positive integer. The index the computed column has in the delay matrix
#' @param N positive integer. Size of delay matrix.
#' @param return_fitted_distribution boolean. If TRUE, the function also returns the gamma distribution that was fitted to the respective column.
#'
#' @return If \code{return_fitted_distribution = FALSE},
#' returns the \code{col_number}th column of the delay matrix, based on the vector of report delays given.
#' If \code{return_fitted_distribution = TRUE}, it returns a list with two elements:
#' \code{column} - delay matrix column as described above,
#' and \code{distribution} - the delay distribution that was fitted to the column.
.get_delay_matrix_column <- function(recent_counts_distribution, fit = "none", col_number, N, return_fitted_distribution = FALSE) {
  .are_valid_argument_values(list(
    list(recent_counts_distribution, "numeric_vector"),
    list(fit, "delay_matrix_column_fit"),
    list(col_number, "positive_integer"),
    list(N, "positive_integer"),
    list(return_fitted_distribution, "boolean")
  ))
  i <- col_number
  new_column <- c()

  if (fit == "gamma") {
    gamma_fit <- try(suppressWarnings(fitdistrplus::fitdist(recent_counts_distribution + 1, distr = "gamma")), silent = T)
    if ("try-error" %in% class(gamma_fit)) {
      # TODO only output this if verbose output
      cat("    mle failed to estimate the parameters. Trying method = \"mme\"\n")
      gamma_fit <- fitdistrplus::fitdist(recent_counts_distribution + 1, distr = "gamma", method = "mme")
      # TODO if none work revert to empirical distribution
    }

    shape_fit <- gamma_fit$estimate[["shape"]]
    scale_fit <- 1 / gamma_fit$estimate[["rate"]]

    distribution <- list(name = "gamma", shape = shape_fit, scale = scale_fit)
    delay_distr <- build_delay_distribution(distribution, offset_by_one = TRUE)
  } else { # no fit
    delay_distr <- graphics::hist(recent_counts_distribution, breaks = seq(0, N, l = N + 1), plot = FALSE)
    delay_distr <- delay_distr$density
  }

  if (length(delay_distr) < N - i + 1) {
    delay_distr <- c(delay_distr, rep(0, times = N - i + 1 - length(delay_distr)))
  }
  new_column <- c(rep(0, times = i - 1), delay_distr[1:(N - i + 1)])

  if (fit == "gamma" && return_fitted_distribution == TRUE) {
    return(list(column = new_column, distribution = distribution))
  }
  return(new_column)
}


#' Utility function that generates delay data, assuming a different delay between event and observation for each individual day.
#' It then generates the delay matrix and computes the RMSE between the parameters of the gamma distributions passed as arguments and the ones recovered from the delay matrix.
#' The shapes and scales of the gamma distributions are specified as parameters, and the number of timesteps is assumed to be equal to the length of these vectors.
#'
#' @param original_distribution_shapes vector. Specifies the shapes for the gamma distributions.
#' @param original_distribution_scales vector. Specifies the scales for the gamma distributions.
#' @param nr_distribution_samples integer. How many cases to be sampled for each timestep.
#'
#' @return A list with the computed RMSE. It has two elements: $shape_rmse and $scale_rmse
.delay_distribution_matrix_rmse_compute <- function(original_distribution_shapes, original_distribution_scales, nr_distribution_samples = 500) {

  # Create a vector with all dates in observation interval
  start_date <- as.Date("2021/04/01")
  time_steps <- length(original_distribution_shapes)
  end_date <- start_date + time_steps
  available_dates <- seq(start_date, end_date, by = "day")

  # Build the delay data; Events on each individual day are assumed to be observed according to a different gamma distribution, as specified by original_distribution_shapes and original_distribution_scales,
  sampled_report_delays <- c()
  report_dates <- as.Date(c())
  for (i in 1:time_steps) {
    new_sampled_report_delays <- .sample_from_distribution(list(name = "gamma", shape = original_distribution_shapes[i], scale = original_distribution_scales[i]), nr_distribution_samples)
    sampled_report_delays <- c(sampled_report_delays, new_sampled_report_delays)
    new_report_dates <- rep(available_dates[i], nr_distribution_samples)
    report_dates <- c(report_dates, new_report_dates)
  }
  delay_data <- dplyr::tibble(event_date = report_dates, report_delay = sampled_report_delays)
  result <- get_matrix_from_empirical_delay_distr(delay_data, time_steps, fit = "gamma", return_fitted_distribution = TRUE)

  delay_matrix <- result$matrix
  distrib_list <- result$distributions


  # Get the shapes and scales of the gamma distributions fitted by the get_matrix_from_empirical_delay_distr function
  distribution_shapes <- c()
  distribution_scales <- c()

  for (distribution in distrib_list) {
    distribution_shapes <- c(distribution_shapes, distribution$shape)
    distribution_scales <- c(distribution_scales, distribution$scale)
  }

  # Compute the RMSE between the desired gamma distribution shapes and scales, and the ones obtained by the get_matrix_from_empirical_delay_distr function
  start_index <- length(distribution_shapes) - length(original_distribution_shapes) + 1
  shape_rmse <- Metrics::rmse(distribution_shapes[start_index:length(distribution_shapes)], original_distribution_shapes) / mean(original_distribution_shapes)
  scale_rmse <- Metrics::rmse(distribution_scales[start_index:length(distribution_scales)], original_distribution_scales) / mean(original_distribution_scales)

  return(list(shape_rmse = shape_rmse, scale_rmse = scale_rmse))
}
