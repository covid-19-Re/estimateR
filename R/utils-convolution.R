
#TODO add input validations in all utilities (create validation utilities when needed)


#TODO fill documentation
#' Title
#'
#' @param vector_a numeric
#' @param vector_b numeric
#'
#' @return numeric vector
.convolve_delay_distribution_vectors <- function(vector_a, vector_b){

  .check_is_probability_distr_vector(vector_a)
  .check_is_probability_distr_vector(vector_b)

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

  #TODO add check that matrix_b is lower-triangular (add utility for that)

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

#TODO generalize to any number of delays
#TODO all for list of distributions as input
#TODO ugly if-else structure: improve if possible
#TODO fill doc
#' Title
#'
#' @param first_delay
#' @param second_delay
#'
#' @return vector or matrix
.convolve_delay_distributions <- function(first_delay,
                                          second_delay) {

  #TODO check that elements are either delay distribution vectors or matrices

  if( identical( class(first_delay), "numeric") ){
    if (identical( class(second_delay), "numeric")){
      return(.convolve_delay_distribution_vectors(first_delay, second_delay))
    } else if ( "matrix" %in% class(second_delay) ) {
      return(.convolve_delay_distribution_vector_with_matrix(first_delay, second_delay, vector_first = TRUE))
    } else {
      stop("'second_delay' must be a numeric vector or matrix.")
    }
  } else if( "matrix" %in% class(first_delay) ) {
    if ( identical( class(second_delay), "numeric") ){
      return(.convolve_delay_distribution_vector_with_matrix(first_delay, second_delay, vector_first = FALSE))
    } else if ( "matrix" %in% class(second_delay) ) {
      #TODO work on .convolve_delay_distribution_matrices()
      stop("Convolution of two matrices is not properly implemented")
      # return(.convolve_delay_distribution_matrices(first_delay, second_delay))
    } else {
      stop("'second_delay' must be a numeric vector or matrix.")
    }
  } else {
    stop("'first_delay' must be a numeric vector or matrix.")
  }
}


#TODO generalize this to a list of inputs
#TODO fill in doc
#' Title
#'
#' @param delay_incubation
#' @param delay_onset_to_report
#' @param n_report_time_steps Integer value required if
#' @param start_date
#' @param time_step
#' @param min_number_cases
#' @param upper_quantile_threshold
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
  } else {
    delay_distribution_incubation <- .get_delay_distribution(delay_incubation)
  }

  if( .check_is_empirical_delay_data(delay_onset_to_report) ){
    if(n_report_time_steps == 0) {
      stop("Empirical delay data input but 'n_time_steps' parameter was not set or set to zero.")
    }
    delay_distribution_onset_to_report <- get_matrix_from_empirical_delay_distr(empirical_delays = delay_onset_to_report,
                                                                                n_report_time_steps = n_report_time_steps,
                                                                                start_date = start_date,
                                                                                time_step = time_step,
                                                                                min_number_cases = min_number_cases,
                                                                                upper_quantile_threshold = upper_quantile_threshold)
  } else {
    delay_distribution_onset_to_report <- .get_delay_distribution(delay_onset_to_report)
  }

  total_delay_distribution <- .convolve_delay_distributions(delay_distribution_incubation,
                                                            delay_distribution_onset_to_report)

  return(total_delay_distribution)
}
