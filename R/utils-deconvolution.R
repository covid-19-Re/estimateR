#TODO write deconvolution step with the convolution of a matrix with vector. That should apply to Swiss non-onset data
#TODO write deconvolution step that takes into account that data is onset data (when it is). That should apply to Swiss onset data

#TODO add a function to deal with empirical delay distribution and build delay distribution matrix
## make sure it can deal with Spanish data specificity

#TODO transform make_ecdf_from_two_gammas into a more general function summing draws from n distributions (allow type of distribution to be changed)

#TODO make make_ecdf_from_empirical_data_and_gamma more general, draws can be from different types of distribution.
## Check why gamma_draws is an input and not drawn inside the function
## Check why not Vectorized ecdf output as in make_ecdf_from_two_gammas

#TODO think about whether utilities like '.get_matrix_from_single_delay_distr' need to be exported


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


#TODO fill documentation
#' Title
#'
#' @param vector_a numeric
#' @param vector_b numeric
#'
#' @return numeric vector
.convolve_delay_distribution_vectors <- function(vector_a, vector_b){
  # Right-pad vectors with zeroes to bring them to the same length

  final_length <- length(vector_b) + length(vector_a)

  vector_a <- c(vector_a, rep(0,times = final_length - length(vector_a)))
  vector_b <- c(vector_b, rep(0,times = final_length - length(vector_b)))

  vector_c <- rep(0, times = final_length)

  for(i in 1 : final_length) {
    reversed_vector_b <- rev(vector_b[1:i]) # Reverse vector_b truncated at index i
    reversed_vector_b <- c(reversed_vector_b, rep(0, times = final_length - i)) # Right-pad with zeroes
    vector_c[i] <- vector_a %*% reversed_vector_b # Compute dot product between vectors
  }

  return(vector_c)
}

# TODO document function and document delay distribution matrix format
# TODO if vector_first = FALSE is used, need to consider doing an operation equivalent to .left_augment_delay_distribution
#' Title
#'
#' @param vector_a
#' @param matrix_b
#' @param vector_first a boolean. Delay described in vector is applied before delay described in matrix
#'
#' @return square matrix. Delay distribution matrix
.convolve_delay_distribution_vector_with_matrix <- function(vector_a, matrix_b, vector_first = TRUE){

  #TODO add check that matrix_b is square

  if( vector_first ) {
    n_col_augment <- .get_initial_deconvolution_shift(vector_a)
    matrix_b <- .left_augment_delay_distribution(delay_distribution_matrix = matrix_b,
                                                 n_col_augment = n_col_augment)
  }

   N <- nrow(matrix_b)
   # Right-pad vector with zeroes to bring to same dimension as square matrix
   vector_a <- c(vector_a, rep(0, times = max(0, N - length(vector_a))))

   # Initialize result matrix
   convolved_matrix <- matrix(0, nrow = N, ncol = N)

   # Iterate over columns (each column represents the delay distribution on a specific date)
   for(j in 1:N) {
     # Iterate over rows
      for(i in 0 : (N - j)) {
        if(vector_first) { # Take corresponding row in matrix_b
          # The row is left-truncated (only j to N indices) so as to start at same date (date with index j) as column in convolved matrix
          matrix_b_elements <- matrix_b[i + j, j : (j + i) ]
        } else { # Take corresponding column in matrix_b (and revert it)
          matrix_b_elements <- matrix_b[(i + j) : j, j]
        }

        truncated_vector_a <- vector_a[1:(i+1)]
        convolved_matrix[i + j, j] <- truncated_vector_a %*% matrix_b_elements
      }
   }

   return(convolved_matrix)
}

# TODO document function and document delay distribution matrix format
# TODO if used, need to consider if left augmentation is required like for .convolve_delay_distribution_vector_with_matrix
#' Title
#'
#' Note that this convolution operation is not commutative!
#' @param matrix_a square numeric matrix
#' @param matrix_b square numeric matrix
#'
#' @return square matrix. Convolved matrix of time-varying delays.
.convolve_delay_distribution_matrices <- function(matrix_a, matrix_b){
  #TODO return error if matrices are not square or not of the same size

  N <- nrow(matrix_a)
  # Initialize result matrix
  convolved_matrix <- matrix(0, nrow = N, ncol = N)

  # Iterate over columns (each column represents the delay distribution on a specific date)
  for(j in 1:N) {
    # Iterate over rows
    for(i in 0 : (N - j)) {

      # Take truncated column of matrix_a (first delay applied)
      matrix_a_elements <- matrix_a[ j : (j + i), j ]
      # Take truncated row of matrix_b (second delay applied)
      matrix_b_elements <- matrix_b[i + j, j : (j + i) ]

      convolved_matrix[i + j, j] <- matrix_a_elements %*% matrix_b_elements
    }
  }
  return(convolved_matrix)

}


#TODO fill in and update documentation
#' Build a waiting time distribution from the convolution of two gamma distributions
#'
#' @param distribution_incubation list. delay between infection and symptom onset
#' @param distribution_onset_to_report list. delay between symtom onset and observation
#' @param max_quantile numeric. between 0 and 1.
#'
#' @return vector specifying the CDF between each time step of the waiting time distribution.
#' @export
#'
#' @examples
#' #TODO add examples
combine_incubation_with_reporting_delay <- function(distribution_incubation,
                                                    distribution_onset_to_report,
                                                    max_quantile = 0.9999,
                                                    incubation_offset_by_one = FALSE,
                                                    onset_to_report_offset_by_one = FALSE) {


  delay_distribution_incubation <- build_delay_distribution(distribution = distribution_incubation,
                                                            max_quantile = max_quantile,
                                                            offset_by_one = incubation_offset_by_one)

  delay_distribution_onset_to_report <- build_delay_distribution(distribution = distribution_onset_to_report,
                                                                 max_quantile = max_quantile,
                                                                 offset_by_one = onset_to_report_offset_by_one)


  convolved_output <- .convolve_delay_distribution_vectors(delay_distribution_incubation,
                                                           delay_distribution_onset_to_report)

  return(convolved_output)
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
#' @param time_step string. "day", "X days", "week", "month"... (see \link{\code{seq.Date}} for details)
#' @param min_number_cases integer. Minimal number of cases to build empirical distribution from
#' @param upper_quantile_threshold numeric. Between 0 and 1. TODO add details
#'
#' @return
#' @export
#'
#' @examples
#' #TODO add example
get_matrix_from_empirical_delay_distr <- function(empirical_delays,
                                                  start_date,
                                                  n_report_time_steps,
                                                  time_step = "day",
                                                  min_number_cases = 300,
                                                  upper_quantile_threshold = 0.99){

  ##TODO need to account for offset if onset data (or not onset data?)
  ##TODO reconsider if we make gamma fit (allow to turn it off, or to use different distribution)

  all_report_dates <- seq.Date(from = start_date, by = time_step, length.out = n_report_time_steps)

  # Ignore the delay data that is posterior to the last incidence report date.
  empirical_delays <- empirical_delays %>%
    dplyr::filter(event_date <= max(all_report_dates))

  # Find the threshold for right-truncation
  # No time-variation beyond this threshold due to the fraction of unsampled individuals when nearing the last sampling date
  # TODO put the search for threshold_right_truncation in separate utility function
  delay_counts <- empirical_delays %>%
    dplyr::select(report_delay) %>%
    dplyr::group_by(report_delay) %>%
    dplyr::summarise(counts = dplyr::n(), .groups = "drop")

  threshold_right_truncation <- delay_counts %>%
    dplyr::mutate(cumul_freq = cumsum(counts)/sum(counts)) %>%
    dplyr::filter(cumul_freq > upper_quantile_threshold) %>%
    utils::head(n=1) %>%
    dplyr::pull(report_delay)

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
      dplyr::arrange( dplyr::desc(event_date) ) %>%
      dplyr::filter( event_date <= all_dates[i] )

    if( nrow(recent_counts) >= min_number_cases ) {
      # If enough data points before date of interest,
      # take most recent observations before this date

      recent_counts_distribution <- recent_counts %>%
        dplyr::slice_head( n = min_number_cases )  %>%
        dplyr::pull(report_delay)
    } else {
      # Otherwise, take 'min_number_of_cases' observations,
      # even after date of interest.
      recent_counts_distribution <- shuffled_delays %>%
        dplyr::arrange( event_date ) %>%
        dplyr::slice_head( n = min_number_cases )  %>%
        dplyr::pull(report_delay)
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


