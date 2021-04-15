#TODO reorganize files between utils-deconvolution, utils-convolution, utils-distribution

#TODO write deconvolution step with the convolution of a matrix with vector. That should apply to Swiss non-onset data

#TODO write deconvolution step that takes into account that data is onset data (when it is). That should apply to Swiss onset data

#TODO add a way to deal with Spanish data specificity

#TODO add way to pass list of distributions in pipe

#' Make square delay distribution matrix from vector of delay distributions.
#'
#' @param delay_distribution numeric vector
#' @param N integer. Dimension of output matrix
#'
#' @return square numeric matrix
.get_matrix_from_single_delay_distr <- function(delay_distribution, N) {

  if(N >= length(delay_distribution)) {
    delay_distribution <- c(delay_distribution, rep(0, times = N - length(delay_distribution)))
  }

  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), delay_distribution[1:(N - i + 1)])
  }

  return(delay_distribution_matrix)
}

#TODO fill documentation
#TODO maybe export
#TODO maybe merge with .get_matrix_from_single_delay_distr by adding N parm and checking if list or unique vector
#TODO finish replacing build_delay_distribution input.
#' Build delay distribution matrix from list of delay distributions
#'
#' @param distributions list of distributions
#' @param max_quantile
#' @param offset_by_one A boolean
#'
#' @return delay distribution matrix
.get_delay_matrix_from_delay_distribution_parms <- function(distributions,
                                                            max_quantile = 0.999,
                                                            offset_by_one = FALSE) {
  #TODO add checks on validity of input


  # Generate list of delay distribution vectors
  delay_distribution_list <- lapply(distributions, function(distr){
    build_delay_distribution(distr,
                            max_quantile = max_quantile,
                            offset_by_one = offset_by_one)
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

#TODO document
#' Title
#'
#' @param delay_distribution_matrix
#' @param n_col_augment an integer.
#'
#' @return
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

#TODO document
#TODO add options to take median, mode, mean...
#' Get initial shift for deconvolution step
#'
#' @param delay_distribution_vector
#'
#' @return
.get_initial_deconvolution_shift <- function(delay_distribution_vector){
    initial_shift <- ceiling(min(which(cumsum(delay_distribution_vector) > 0.5)) - 1)
    return(initial_shift)
}




#TODO test
#TODO improve function documentation
#TODO possibly allow for other ways to specifiy initial sift than median of all reports.
#TODO format of empirical_delays must be specified somewhere:
# use "event_date" and "report_delay" as column names
#TODO allow for other distributions than gamma for fit, and also allow no fit.

#' Build matrix of delay distributions through time from empirical delay data.
#'
#' This matrix is required for the application of the Richardson-Lucy algorithm.
#'
#' @param empirical_delays tibble. format to be specified
#' @param start_date Date. First date of incidence data
#' @param n_report_time_steps integer. Length of incidence time series
#' @param time_step string. "day", "X days", "week", "month"... (see \code{\link[base]{seq.Date}} for details)
#' @param min_number_cases integer. Minimal number of cases to build empirical distribution from
#' @param upper_quantile_threshold numeric. Between 0 and 1. TODO add details
#'
#' @return
#' @export
#'
#' @examples
#' #TODO add example
get_matrix_from_empirical_delay_distr <- function(empirical_delays,
                                                  n_report_time_steps,
                                                  start_date = NULL,
                                                  time_step = "day",
                                                  min_number_cases = NULL,
                                                  upper_quantile_threshold = 0.99,
                                                  min_number_cases_fraction = 0.05,
                                                  min_min_number_cases = 10){

  ##TODO need to account for offset if onset data (or not onset data?)
  ##TODO reconsider if we make gamma fit (allow to turn it off, or to use different distribution)

  if(is.null(start_date)) {
    start_date <- min(dplyr::pull(empirical_delays, .data$event_date), na.rm = TRUE)
  }

  all_report_dates <- seq.Date(from = start_date, by = time_step, length.out = n_report_time_steps)

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

  # Use median of reported delays as initial shift (needed for deconvolution step)
  initial_shift <- round(stats::median(empirical_delays$report_delay, na.rm = T))

  # Left-pad the dates we are looking at to account for shift between event dates and observation dates.
  all_dates <- c(rev(seq.Date(from = start_date, by = paste0("-1 ", time_step), length.out = initial_shift + 1)),
                 seq.Date(from = start_date, by = time_step, length.out = n_report_time_steps)[-1])

  n_time_steps <- n_report_time_steps + initial_shift

  delay_distribution_matrix <- matrix(0, nrow = n_time_steps, ncol = n_time_steps)

  #TODO fix issue with what happens when n_time_steps <= threshold_right_truncation: we shouldn't go until n_time_steps
  #TODO test what happens when n_time_steps <= threshold_right_truncation
  last_varying_col <- ifelse(n_time_steps > threshold_right_truncation, n_time_steps - threshold_right_truncation, n_time_steps)

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

