
#TODO add input validations in all utilities (create validation utilities when needed)


#' Convolve two discretized probability distribution vectors.
#'
#' @param vector_a,vector_b Discretized probability distribution vectors
#'
#' @return discretized probability distribution vector
.convolve_delay_distribution_vectors <- function(vector_a, vector_b){

  .check_is_probability_distr_vector(vector_a)
  .check_is_probability_distr_vector(vector_b)

  # Right-pad vectors with zeroes to bring them to the same length
  final_length <- length(vector_b) + length(vector_a)
  vector_a <- c(vector_a, rep(0,times = final_length - length(vector_a)))
  vector_b <- c(vector_b, rep(0,times = final_length - length(vector_b)))

  # Initialize result vector
  vector_c <- rep(0, times = final_length)

  # Fill result vector
  for(i in 1 : final_length) {
    reversed_vector_b <- rev(vector_b[1:i]) # Reverse vector_b truncated at index i
    reversed_vector_b <- c(reversed_vector_b, rep(0, times = final_length - i)) # Right-pad with zeroes
    vector_c[i] <- vector_a %*% reversed_vector_b # Compute dot product between vectors
  }

  return(vector_c)
}

# TODO document delay distribution matrix format
# TODO if vector_first = FALSE is used, need to consider doing an operation equivalent to .left_augment_delay_distribution
#' Convolve a delay distribution vector with a delay distribution matrix
#'
#' @param vector_a discretized delay distribution vector
#' @param matrix_b discretized delay distribution matrix
#' @param vector_first Is the delay described by \code{vector_a} happen before
#'   the one from \code{matrix_b} or is it the opposite?
#'
#' @return discretized delay distribution matrix
.convolve_delay_distribution_vector_with_matrix <- function(vector_a, matrix_b, vector_first = TRUE){

  .check_is_probability_distr_vector(vector_a)
  #TODO test that matrix_b is a delay distribution matrix

  if( vector_first ) {
    # Increase size of matrix_b to account for the fact that the output matrix will be shifted in time by the vector_a delay
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
        # The row is left-truncated (only j to N indices) so as to start at same
        # date (date with index j) as column in convolved matrix
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

# TODO if used, need to consider if left augmentation is required like for .convolve_delay_distribution_vector_with_matrix
#' Convolve two delay distribution matrices
#'
#' Note that this convolution operation is not commutative!
#' The order matters: here the delay implied by matrix_a happens first.
#' @param matrix_a first delay distribution matrix
#' @param matrix_b second delay distribution matrix
#'
#' @return convolved discretized delay distribution matrix
.convolve_delay_distribution_matrices <- function(matrix_a, matrix_b){

  stop("This function is not ready.\n
       Need to consider if left augmentation is required like for .convolve_delay_distribution_vector_with_matrix")

  if(nrow(matrix_a) != nrow(matrix_b)) {
    stop("Convolved matrices must have the same dimensions.")
  }

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

#TODO generalize to any number of delays
#TODO ugly if-else structure: improve if possible
#' Convolve two delay vectors or matrices
#'
#' @param first_delay discretized delay distribution vector or matrix
#' @param second_delay discretized delay distribution vector or matrix
#'
#' @return discretized delay distribution vector (if both input delays are
#'   vectors) or matrix
.convolve_delay_distributions <- function(first_delay,
                                          second_delay) {

  #TODO check that elements are either delay distribution vectors or matrices

  if( .is_numeric_vector(first_delay) ){
    if ( .is_numeric_vector(second_delay) ){
      return(.convolve_delay_distribution_vectors(first_delay, second_delay))
    } else if ( is.matrix(second_delay) ) {
      return(.convolve_delay_distribution_vector_with_matrix(vector_a = first_delay,
                                                             matrix_b = second_delay,
                                                             vector_first = TRUE))
    } else {
      stop("'second_delay' must be a numeric vector or matrix.")
    }
  } else if( is.matrix(first_delay) ) {
    if ( .is_numeric_vector(second_delay) ){
      return(.convolve_delay_distribution_vector_with_matrix(matrix_b = first_delay,
                                                             vector_a = second_delay,
                                                             vector_first = FALSE))
    } else if ( is.matrix(second_delay) ) {
      #TODO work on .convolve_delay_distribution_matrices()
      stop("Convolution of two matrices is not properly implemented yet.")
      # return(.convolve_delay_distribution_matrices(first_delay, second_delay))
    } else {
      stop("'second_delay' must be a numeric vector or matrix.")
    }
  } else {
    stop("'first_delay' must be a numeric vector or matrix.")
  }
}

#TODO allow for list of distributions as input
#TODO generalize this to a list of inputs
#TODO fill in doc
#' Title
#'
#' @param delay_incubation Delay
#' @param delay_onset_to_report
#' @param n_report_time_steps Integer value required if
#' @param start_date Optional. Date of the first incidence data point.
#' @param time_step Size of time step used in delays. String in the format:
#'   "day", "2 days", "week", "year"... (see \code{\link[base]{seq.Date}} for
#'   details)
#' @param min_number_cases Use only if an input is empirical delay data. Number
#'   of cases over which to build empirical delay distribution
#' @param upper_quantile_threshold Value between 0 and 1. Quantile
#'
#' @return
#' @export
#'
#' @examples
#' #TODO add example
convolve_delay_inputs <- function(delay_incubation,
                                  delay_onset_to_report,
                                  n_report_time_steps = 0,
                                  start_date = NULL,
                                  time_step = "day",
                                  min_number_cases = NULL,
                                  upper_quantile_threshold = 0.99){

  #TODO put these tests below in a utility function
  if( .check_is_empirical_delay_data(delay_incubation) ){
    if(n_report_time_steps == 0) {
      stop("Empirical delay data input but 'n_time_steps' parameter was not set or set to zero.")
    }
    delay_distribution_incubation <- get_matrix_from_empirical_delay_distr(empirical_delays = delay_incubation,
                                                                           n_report_time_steps = n_report_time_steps,
                                                                           start_date = start_date,
                                                                           time_step = time_step,
                                                                           min_number_cases = min_number_cases,
                                                                           upper_quantile_threshold = upper_quantile_threshold)
  } else if(is.matrix(delay_incubation)) {
    #TODO check that delay_incubation is actually a delay distribution matrix
    delay_distribution_incubation <- delay_incubation
  } else if(.is_numeric_vector(delay_incubation) || is.list(delay_incubation)){
    delay_distribution_incubation <- .get_delay_distribution(delay_incubation)
  } else {
    stop("Invalid 'delay_incubation' input. 'delay_incubation' must be either:
         a numeric vector representing a discretized probability distribution,
         or a matrix representing discretized probability distributions,
         or a distribution object (e.g. list(name = 'gamma', scale = X, shape = Y)),
         or empirical delay data.")
  }

  if( .check_is_empirical_delay_data(delay_onset_to_report) ){
    if(n_report_time_steps == 0) {
      stop("Empirical delay data input but 'n_report_time_steps' parameter was not set or set to zero.")
    }
    delay_distribution_onset_to_report <- get_matrix_from_empirical_delay_distr(empirical_delays = delay_onset_to_report,
                                                                                n_report_time_steps = n_report_time_steps,
                                                                                start_date = start_date,
                                                                                time_step = time_step,
                                                                                min_number_cases = min_number_cases,
                                                                                upper_quantile_threshold = upper_quantile_threshold)
  } else if(is.matrix(delay_onset_to_report)) {
    #TODO check that delay_incubation is actually a delay distribution matrix
    delay_distribution_onset_to_report <- delay_onset_to_report
  } else if (.is_numeric_vector(delay_onset_to_report) || is.list(delay_onset_to_report)) {
    delay_distribution_onset_to_report <- .get_delay_distribution(delay_onset_to_report)
  } else {
    stop("Invalid 'delay_onset_to_report' input. 'delay_onset_to_report' must be either:
         a numeric vector representing a discretized probability distribution,
         or a matrix representing discretized probability distributions,
         or a distribution object (e.g. list(name = 'gamma', scale = X, shape = Y)),
         or empirical delay data.")
  }

  total_delay_distribution <- .convolve_delay_distributions(delay_distribution_incubation,
                                                            delay_distribution_onset_to_report)

  return(total_delay_distribution)
}
