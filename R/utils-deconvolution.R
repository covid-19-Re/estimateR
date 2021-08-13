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

  .are_valid_argument_values(list(list(delay_distribution_vector, "probability_distr_vector"),
                                  list(N, "positive_integer")))

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

  for(i in 1:length(distributions)){
    .are_valid_argument_values(list(list(distributions[[i]], "distribution")))
  }
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

  .are_valid_argument_values(list(list(delay_distribution_matrix, "probability_distr_matrix", 0),
                                  list(n_col_augment, "positive_integer")))

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
#TODO validate 'quantile'
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
  .are_valid_argument_values(list(list(delay_distribution_vector, "probability_distr_vector")))

    initial_shift <- ceiling(min(which(cumsum(delay_distribution_vector) > quantile))) - 1
    initial_shift <- max(initial_shift, 0, na.rm = TRUE)
    return(initial_shift)
}


#TODO test
#TODO possibly allow for other ways to specifiy initial shift than median of all reports.
#TODO allow for other distributions than gamma for fit.
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
#' @param fit string. One of "gamma" or "none". Specifies the type of fit that 
#' is applied to the columns of the delay matrix
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
                                                  min_number_cases_fraction = 0.05,
                                                  min_min_number_cases = 10,
                                                  fit = "none",
                                                  return_fitted_distribution = FALSE){
  ##TODO need to account for offset if onset data (or not onset data?)
  ##TODO reconsider if we make gamma fit (allow to turn it off, or to use different distribution)

  .are_valid_argument_values(list(list(empirical_delays, "empirical_delay_data"),
                                  list(n_report_time_steps, "positive_integer"),
                                  list(ref_date, "null_or_date"),
                                  list(time_step, "time_step"),
                                  list(min_number_cases, "null_or_int"),
                                  list(upper_quantile_threshold, "numeric_between_zero_one"),
                                  list(min_number_cases_fraction, "numeric_between_zero_one"),
                                  list(min_min_number_cases, "positive_integer"),
                                  list(fit, "delay_matrix_column_fit"),
                                  list(return_fitted_distribution, "boolean")))

  
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
  #TODO it may be simpler to just do the augmentation during the deconvolution step
  initial_shift <- ceiling(stats::quantile(empirical_delays$report_delay, probs = 0.99, na.rm = T))[1]

  # Left-pad the dates we are looking at to account for shift between event dates and observation dates.
  all_dates <- c(rev(seq.Date(from = ref_date, by = paste0("-1 ", time_step), length.out = initial_shift + 1)),
                 seq.Date(from = ref_date, by = time_step, length.out = n_report_time_steps)[-1])

  n_time_steps <- n_report_time_steps + initial_shift

  delay_distribution_matrix <- matrix(0, nrow = n_time_steps, ncol = n_time_steps)

  last_varying_col <- dplyr::if_else(n_time_steps > threshold_right_truncation, n_time_steps - threshold_right_truncation, n_time_steps)
  
  distrib_list <- list() #needed for the test that checks if get_matrix_from_empirical_delay_distr returns a matrix with the expected distributions when using fit = "gamma"

      # Populate the delay_distribution_matrix by column
    if(n_time_steps > threshold_right_truncation){
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

        result <- .get_delay_matrix_column(recent_counts_distribution, fit, col_number = i, n_time_steps, return_fitted_distribution)
        if(is.list(result)){
          distrib_list[[i]] <- result$distribution
          new_column <- result$column
        } else {
          new_column <- result
        }
        delay_distribution_matrix[,i] <- new_column
      }
    } else { # if n_time_steps <= threshold_right_truncation
        # Shuffle rows so as to get rid of potential biases
        shuffled_delays <- empirical_delays %>%
          dplyr::slice( sample(1:dplyr::n()) )

        recent_counts <- shuffled_delays %>%
          dplyr::arrange( dplyr::desc(.data$event_date) ) %>%
          dplyr::filter( .data$event_date <= all_dates[1] )

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

        result <- .get_delay_matrix_column(recent_counts_distribution, fit, col_number = 1, n_time_steps, return_fitted_distribution)
        if(is.list(result)){
          distrib_list[[i]] <- result$distribution
          new_column <- result$column
        } else {
          new_column <- result
        }
        
        for(i in 0:(last_varying_col-1)) {
          delay_distribution_matrix[, i+1] <- c(rep(0, times =  i), new_column[1:(length(new_column) - i)])
          if(fit == "gamma"){
            distrib_list <- append(distrib_list, distrib_list[length(distrib_list)])
          }
        }
    }
  
    if(last_varying_col < n_time_steps) {
      for( j in 1: threshold_right_truncation ) {
          delay_distribution_matrix[, i+j ] <-  c(rep(0, times =  j), delay_distribution_matrix[1:(nrow(delay_distribution_matrix) - j), i])
          if(fit == "gamma"){
            distrib_list <- append(distrib_list, distrib_list[length(distrib_list)])
          }
      }
    }
  
  if(return_fitted_distribution){
    return(list(matrix = delay_distribution_matrix, distributions = distrib_list))
  }
  return( delay_distribution_matrix )
}


#' Build a specific column of the delay distribution matrix
#'
#' @param recent_counts_distribution numeric vector of report delays, as used in \code{get_matrix_from_empirical_delay_distr}
#' @param fit string. Can be either "none" or "gamma". Specifies the type of fitting applied to the computed column
#' @param col_number positive integer. The index the computed column has in the delay matrix
#' @param N positive integer. Size of delay matrix.
#' @param return_fitted_distribution boolean. If TRUE, the function also returns the gamma distribution that was fitted to the respective column.
#'
#' @return If \code{return_fitted_distribution = FALSE}, returns the \code{col_number}th column of the delay matrix, based on the vector of report delays given. 
#' If \code{return_fitted_distribution = TRUE}, it returns a list with two elements: \code{column} - delay matrix column as described above, and \code{distribution} - the delay distribution that was fitted to the column. 
.get_delay_matrix_column <- function(recent_counts_distribution, fit = "none", col_number, N, return_fitted_distribution = FALSE){
  .are_valid_argument_values(list(list(recent_counts_distribution, "numeric_vector"),
                                  list(fit, "delay_matrix_column_fit"),
                                  list(col_number, "positive_integer"),
                                  list(N, "positive_integer"),
                                  list(return_fitted_distribution, "boolean")))
  i <- col_number
  new_column <- c()

  if(fit == "gamma"){
    gamma_fit <- try(suppressWarnings(fitdistrplus::fitdist(recent_counts_distribution + 1, distr = "gamma")), silent = T)
    if ("try-error" %in% class(gamma_fit)) {
      #TODO only output this if verbose output
      cat("    mle failed to estimate the parameters. Trying method = \"mme\"\n")
      gamma_fit <- fitdistrplus::fitdist(recent_counts_distribution + 1, distr = "gamma", method = "mme")
    }
    #TODO if none work revert to empirical distribution

    shape_fit <- gamma_fit$estimate[["shape"]]
    scale_fit <- 1/gamma_fit$estimate[["rate"]]

    distribution <- list(name = "gamma", shape = shape_fit, scale = scale_fit)
    delay_distr <- build_delay_distribution(distribution, offset_by_one=TRUE)
    
  } else { #no fit
    delay_distr <- hist(recent_counts_distribution, breaks=seq(0, N, l = N + 1), plot=FALSE)
    delay_distr <- delay_distr$density
  }

  if(length(delay_distr) < N - i + 1) {
    delay_distr <- c(delay_distr, rep(0, times = N - i + 1 - length(delay_distr)))
  }
  new_column <-  c(rep(0, times = i - 1 ), delay_distr[1:(N - i + 1)])
  
  if(fit == "gamma" && return_fitted_distribution == TRUE){
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
.delay_distribution_matrix_rmse_compute <- function(original_distribution_shapes, original_distribution_scales, nr_distribution_samples = 500){

  #Create a vector with all dates in observation interval
  start_date <- as.Date('2021/04/01')
  time_steps = length(original_distribution_shapes)
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



