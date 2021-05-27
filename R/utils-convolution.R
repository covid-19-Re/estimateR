
#TODO add input validations in all utilities (create validation utilities when needed)

#' Convolve two discretized probability distribution vectors.
#'
#' @inheritParams distribution
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

#' Convolve a delay distribution vector with a delay distribution matrix
#'
#' @inheritParams distribution
#' @param vector_first Does the delay described by \code{vector_a} happen before
#'   the one from \code{matrix_b}?
#'
#' @return discretized delay distribution matrix
.convolve_delay_distribution_vector_with_matrix <- function(vector_a, matrix_b, vector_first = TRUE){

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
#' The order matters: here the delay implied by \code{matrix_a} happens first,
#' then the one implied by \code{matrix_b}.
#' @inheritParams distribution
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

  #TODO check that elements are either delay distribution vectors or matrices: validate input

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

#TODO generalize this to a list of inputs
#TODO add details on what is returned (vector or matrix)
#TODO rename to 'convolve_delays'
#' Convolve delay inputs.
#'
#' This function is flexible in the type of delay inputs it can handle.
#' Each delay input can be one of:
#' \itemize{
#' \item{a list representing a distribution object}
#' \item{a discretized delay distribution vector}
#' \item{a discretized delay distribution matrix}
#' \item{a dataframe containing empirical delay data}
#' }
#'
#' see \code{\link{get_matrix_from_empirical_delay_distr}} for details on the format
#' expected for the empirical delay data.
#'
#'
#'
#' @param delay_incubation Incubation delay. Flexible format.
#' @param delay_onset_to_report Delay between symptom onset and report. Flexible format.
#' @param n_report_time_steps integer. Length of incidence time series.
#' Use only when providing empirical delay data.
#' @inheritParams dating
#' @inheritDotParams get_matrix_from_empirical_delay_distr -empirical_delays
#'
#' @return a discretized delay distribution vector or matrix.
#' @export
convolve_delay_inputs <- function(delay_incubation,
                                  delay_onset_to_report,
                                  n_report_time_steps = NULL,
                                  ...){

  .are_valid_argument_values(list(list(delay_incubation, "delay_object", 1),
                                  list(delay_onset_to_report, "delay_object", 1)))

  dots_args <- .get_dots_as_list(...)

  delay_distribution_incubation <- do.call(
    '.get_delay_distribution',
    c(list(delay = delay_incubation,
           n_report_time_steps = n_report_time_steps),
      .get_shared_args(list(get_matrix_from_empirical_delay_distr,
                            build_delay_distribution), dots_args))
  )

  delay_distribution_onset_to_report <- do.call(
    '.get_delay_distribution',
    c(list(delay = delay_onset_to_report,
           n_report_time_steps = n_report_time_steps),
      .get_shared_args(list(get_matrix_from_empirical_delay_distr,
                            build_delay_distribution), dots_args))
  )

  total_delay_distribution <- .convolve_delay_distributions(delay_distribution_incubation,
                                                            delay_distribution_onset_to_report)
  return(total_delay_distribution)
}
