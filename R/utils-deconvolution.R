#TODO reorganize files between utils-deconvolution, utils-convolution, utils-distribution
#TODO write deconvolution step that takes into account that data is onset data (when it is). That should apply to Swiss onset data
#TODO add a way to deal with Spanish data specificity
#TODO add way to pass list of distributions in pipe

#' Make delay distribution matrix from vector of delay distribution.
#'
#' @inheritParams distribution
#' @param N integer. Dimension of output matrix.
#'
#' @return discretized delay distribution matrix, representing a constant-through-time
#' delay distribution.
.get_matrix_from_single_delay_distr <- function(delay_distribution_vector, N) {

  if(N >= length(delay_distribution_vector)) {
    delay_distribution_vector <- c(delay_distribution_vector, rep(0, times = N - length(delay_distribution_vector)))
  }

  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), delay_distribution_vector[1:(N - i + 1)])
  }

  return(delay_distribution_matrix)
}

#TODO maybe export
#TODO maybe merge with .get_matrix_from_single_delay_distr by adding N parm and checking if list or unique vector
#' Build delay distribution matrix from list of delay distributions
#'
#' @param distributions list of distributions,
#' each element is either a distribution list or discretized probability distribution vector.
#' @inheritDotParams build_delay_distribution -distribution
#'
#' @return delay distribution matrix
.get_delay_matrix_from_delay_distribution_parms <- function(distributions, ...) {
  #TODO add checks on validity of input

  dots_args <- .get_dots_as_list(...)

  # Generate list of delay distribution vectors
  delay_distribution_list <- lapply(distributions, function(distr){
    do.call(
      'build_delay_distribution',
      c(list(distribution = distr),
        .get_shared_args(build_delay_distribution, dots_args))
    )
  })

  N <- length(distributions)

  # Initialize empty matrix
  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)

  # Fill matrix by column
  for(i in 1:N){
    delay_distr <- delay_distribution_list[[i]]

    # Right-pad delay_distr vector with zeroes if needed
    if(length(delay_distr) < N - i + 1) {
      delay_distr <- c(delay_distr, rep(0, times = N - i + 1 - length(delay_distr)))
    }
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), delay_distr[1:(N - i + 1)])
  }

  return(delay_distribution_matrix)

}

#' Augment a delay distribution by left padding with new columns.
#'
#' This function reshapes a discretized delay distribution matrix
#' by left-padding it with \code{n_col_augment} columns.
#' Because the output matrix must also be lower-triangular,
#' additional rows are also padded to the top rows.
#' This function allows one to extend further in the past
#' the range of the initial delay distribution matrix.
#' This is useful when convolving that delay distribution matrix
#' with another delay distribution.
#'
#' The columns that are added replicate the left-most column of
#' \code{delay_distribution_matrix}.
#'
#' @inheritParams distribution
#' @param n_col_augment an integer. Number of columns to left-pad
#' \code{delay_distribution_matrix} with.
#'
#' @return If \code{delay_distribution_matrix} is of dimension N,
#' then the result is of dimension N + \code{n_col_augment}.
.left_augment_delay_distribution <- function(delay_distribution_matrix,
                                             n_col_augment){

  n_col_original <- ncol(delay_distribution_matrix)
  n_col_augmented <- n_col_original + n_col_augment

  # Initialize empty matrix
  augmented_matrix <- matrix(0,nrow = n_col_augmented, ncol = n_col_augmented)

  # Fill matrix by column

  # Start by duplicating first column in original matrix into 'n_col_augment' first columns of augmented_matrix
  for(i in 1:n_col_augment){
    augmented_matrix[, i ] <-  c(rep(0, times = i - 1 ), delay_distribution_matrix[, 1], rep(0, times = n_col_augment - i + 1))
  }

  # Then fill with original matrix, adding the required zero on the top rows
  for(i in (n_col_augment + 1):n_col_augmented){
    augmented_matrix[, i ] <-  c(rep(0, times = n_col_augment ), delay_distribution_matrix[, i - n_col_augment])
  }

  return(augmented_matrix)
}

#TODO redoc
#TODO add options to take median, mode, mean...
#TODO redoc quantile parm
#' Get initial shift for deconvolution step
#'
#' This utility function returns the number of timesteps
#' by which the incidence data should be shifted back in the past
#' for the initial step of the Richardson-Lucy deconvolution algorithm.
#'
#' @inheritParams distribution
#'
#' @return an integer value corresponding to the rounded median of
#' the input delay distribution.
.get_time_steps_quantile <- function(delay_distribution_vector, quantile = 0.5){
    initial_shift <- ceiling(min(which(cumsum(delay_distribution_vector) > quantile))) - 1
    initial_shift <- max(initial_shift, 0, na.rm = TRUE)
    return(initial_shift)
}




#TODO test
#TODO possibly allow for other ways to specifiy initial shift than median of all reports.
#TODO allow for other distributions than gamma for fit, and also allow no fit.
#TODO doc fact that if ref_date is not given then ref_date is earliest reported date.
#TODO also doc emphasize fact that ref_date parameter is important here.
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
#' @inherit empirical_delay_data_format
#'
#' @param empirical_delays dataframe containing the empirical data. See Details.
#' @param n_report_time_steps integer. Length of the incidence time series in the accompanying analysis.
#' This argument is needed to determine the dimensions of the output matrix.
#' @param min_number_cases integer. Minimal number of cases to build
#' the empirical distribution from. TODO add details
#' @param min_number_cases_fraction numeric. Between 0 and 1.
#' If \code{min_number_cases} is not provided (kept to \code{NULL}),
#' the number of most-recent cases used to build
#' the instant delay distribution is \code{min_number_cases_fraction}
#' times the total number of reported delays.
#' @param min_min_number_cases numeric. Lower bound
#' for number of cases used to build an instant delay distribution.
#' @param upper_quantile_threshold numeric. Between 0 and 1. TODO add details
#' @inheritParams dating
#'
#' @return a discretized delay distribution matrix.
#' @export
get_matrix_from_empirical_delay_distr <- function(empirical_delays,
                                                  n_report_time_steps,
                                                  ref_date = NULL,
                                                  time_step = "day",
                                                  min_number_cases = NULL,
                                                  upper_quantile_threshold = 0.99,
                                                  min_number_cases_fraction = 0.05,
                                                  min_min_number_cases = 10){

  ##TODO need to account for offset if onset data (or not onset data?)
  ##TODO reconsider if we make gamma fit (allow to turn it off, or to use different distribution)

  if(is.null(ref_date)) {
    ref_date <- min(dplyr::pull(empirical_delays, .data$event_date), na.rm = TRUE)
  }

  all_report_dates <- seq.Date(from = ref_date, by = time_step, length.out = n_report_time_steps)

  # Ignore the delay data that is posterior to the last incidence report date.
  empirical_delays <- empirical_delays %>%
    dplyr::filter(.data$event_date <= max(all_report_dates))

  # Set the 'min_number_cases' parameter if not set by the user
  #TODO make this 'min_number_cases' depend on the length of the time_series.
  if( is.null(min_number_cases) ){
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
    dplyr::mutate(cumul_freq = cumsum(.data$counts)/sum(.data$counts)) %>%
    dplyr::filter(.data$cumul_freq > upper_quantile_threshold) %>%
    utils::head(n=1) %>%
    dplyr::pull(.data$report_delay)

  #TODO allow for different ways of specifying initial shift
  # Use median of reported delays as initial shift (needed for deconvolution step)
  initial_shift <- ceiling(stats::quantile(empirical_delays$report_delay, probs = 0.95, na.rm = T))[1]

  # Left-pad the dates we are looking at to account for shift between event dates and observation dates.
  all_dates <- c(rev(seq.Date(from = ref_date, by = paste0("-1 ", time_step), length.out = initial_shift + 1)),
                 seq.Date(from = ref_date, by = time_step, length.out = n_report_time_steps)[-1])

  n_time_steps <- n_report_time_steps + initial_shift

  delay_distribution_matrix <- matrix(0, nrow = n_time_steps, ncol = n_time_steps)

  #TODO fix issue with what happens when n_time_steps <= threshold_right_truncation: we shouldn't go until n_time_steps
  #TODO test what happens when n_time_steps <= threshold_right_truncation
  last_varying_col <- dplyr::if_else(n_time_steps > threshold_right_truncation, n_time_steps - threshold_right_truncation, n_time_steps)

  distribution_list <- vector(mode = "list", length = last_varying_col)

  # Populate the delay_distribution_matrix by column
  for(i in 1:last_varying_col) {

    # Shuffle rows so as to get rid of potential biases
    shuffled_delays <- empirical_delays %>%
      dplyr::slice( sample(1:dplyr::n()) )

    recent_counts <- shuffled_delays %>%
      dplyr::arrange( dplyr::desc(.data$event_date) ) %>%
      dplyr::filter( .data$event_date <= all_dates[i] )

    if( nrow(recent_counts) >= min_number_cases ) {
      # If enough data points before date of interest,
      # take most recent observations before this date.

      recent_counts_distribution <- recent_counts %>%
        dplyr::slice_head( n = min_number_cases )  %>%
        dplyr::pull(.data$report_delay)
    } else {
      # Otherwise, take 'min_number_of_cases' observations,
      # even after date of interest.
      recent_counts_distribution <- shuffled_delays %>%
        dplyr::arrange( .data$event_date ) %>%
        dplyr::slice_head( n = min_number_cases )  %>%
        dplyr::pull(.data$report_delay)
    }

    gamma_fit <- try(suppressWarnings(fitdistrplus::fitdist(recent_counts_distribution + 1, distr = "gamma")),
                     silent = T)
    if ("try-error" %in% class(gamma_fit)) {
      #TODO only output this if verbose output
      cat("    mle failed to estimate the parameters. Trying method = \"mme\"\n")
      gamma_fit <- fitdistrplus::fitdist(recent_counts_distribution + 1, distr = "gamma", method = "mme")
    }
    #TODO if none work revert to empirical distribution

    shape_fit <- gamma_fit$estimate[["shape"]]
    scale_fit <- 1/gamma_fit$estimate[["rate"]]

    distribution_list[[i]] <- list(name = "gamma", shape = shape_fit, scale = scale_fit)
  }

  if(last_varying_col <= n_time_steps) {
    for( i in 1: threshold_right_truncation ) {
      distribution_list <- append(distribution_list, distribution_list[last_varying_col])
    }
  }

  delay_distribution_matrix <- .get_delay_matrix_from_delay_distribution_parms(distribution_list,
                                                                               offset_by_one = TRUE)

  return( delay_distribution_matrix )
}

