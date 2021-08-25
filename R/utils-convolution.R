#' Convolve two discretized probability distribution vectors.
#'
#' @inheritParams distribution
#'
#' @return discretized probability distribution vector
.convolve_delay_distribution_vectors <- function(vector_a, vector_b) {
  .are_valid_argument_values(list(
    list(vector_a, "probability_distr_vector"),
    list(vector_b, "probability_distr_vector")
  ))

  # Right-pad vectors with zeroes to bring them to the same length
  final_length <- length(vector_b) + length(vector_a)
  vector_a <- c(vector_a, rep(0, times = final_length - length(vector_a)))
  vector_b <- c(vector_b, rep(0, times = final_length - length(vector_b)))

  # Initialize result vector
  vector_c <- rep(0, times = final_length)

  # Fill result vector
  for (i in 1:final_length) {
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
.convolve_delay_distribution_vector_with_matrix <- function(vector_a, matrix_b, vector_first = TRUE) {
  .are_valid_argument_values(list(
    list(vector_a, "probability_distr_vector"),
    list(matrix_b, "probability_distr_matrix", 0),
    list(vector_first, "boolean")
  ))
  if (vector_first) {
    # Increase size of matrix_b to account for the fact that the output matrix will be shifted in time by the vector_a delay
    n_col_augment <- .get_time_steps_quantile(vector_a, quantile = 0.5)
    matrix_b <- .left_augment_delay_distribution(
      delay_distribution_matrix = matrix_b,
      n_col_augment = n_col_augment
    )
  }

  N <- nrow(matrix_b)
  # Right-pad vector with zeroes to bring to same dimension as square matrix
  vector_a <- c(vector_a, rep(0, times = max(0, N - length(vector_a))))

  # Initialize result matrix
  convolved_matrix <- matrix(0, nrow = N, ncol = N)

  # Iterate over columns (each column represents the delay distribution on a specific date)
  for (j in 1:N) {
    # Iterate over rows
    for (i in 0:(N - j)) {
      if (vector_first) { # Take corresponding row in matrix_b
        # The row is left-truncated (only j to N indices) so as to start at same
        # date (date with index j) as column in convolved matrix
        matrix_b_elements <- matrix_b[i + j, j:(j + i)]
      } else { # Take corresponding column in matrix_b (and revert it)
        matrix_b_elements <- matrix_b[(i + j):j, j]
      }

      truncated_vector_a <- vector_a[1:(i + 1)]
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
.convolve_delay_distribution_matrices <- function(matrix_a, matrix_b) {
  stop("This function is not ready.\n
       Need to consider if left augmentation is required like for .convolve_delay_distribution_vector_with_matrix")
  .are_valid_argument_values(list(
    list(matrix_a, "probability_distr_matrix", 0),
    list(matrix_b, "probability_distr_matrix", 0)
  ))

  if (nrow(matrix_a) != nrow(matrix_b)) {
    stop("Convolved matrices must have the same dimensions.")
  }

  N <- nrow(matrix_a)
  # Initialize result matrix
  convolved_matrix <- matrix(0, nrow = N, ncol = N)

  # Iterate over columns (each column represents the delay distribution on a specific date)
  for (j in 1:N) {
    # Iterate over rows
    for (i in 0:(N - j)) {

      # Take truncated column of matrix_a (first delay applied)
      matrix_a_elements <- matrix_a[j:(j + i), j]
      # Take truncated row of matrix_b (second delay applied)
      matrix_b_elements <- matrix_b[i + j, j:(j + i)]

      convolved_matrix[i + j, j] <- matrix_a_elements %*% matrix_b_elements
    }
  }
  return(convolved_matrix)
}

#' Convolve a list of delay vectors or matrices
#'
#' @param delay_distributions list of discretized delay distribution vector or matrix
#'
#' @return discretized delay distribution vector (if all input delays are
#'   vectors) or matrix
.convolve_delay_distributions <- function(delay_distributions) {

  .are_valid_argument_values(list(
    # We put '1' here, because we do not care here about checking the dimension of the matrix.
    list(delay_distributions, "delay_single_or_list", 1)
  ))

  number_of_delays_in_list <- length(delay_distributions)

  if(number_of_delays_in_list == 1) {
    return(delay_distributions[[1]])
  }

  if(number_of_delays_in_list == 2) {
    return(.convolve_two_delay_distributions(delay_distributions[[1]], delay_distributions[[2]]))
  }

  last_delay <- delay_distributions[[number_of_delays_in_list]]

  delay_distributions[number_of_delays_in_list] <- NULL
  convolved_without_last_delay <- .convolve_delay_distributions(delay_distributions)

  return(.convolve_two_delay_distributions(convolved_without_last_delay, last_delay))
}


#' Convolve two delay vectors or matrices
#'
#' @param first_delay discretized delay distribution vector or matrix
#' @param second_delay discretized delay distribution vector or matrix
#'
#' @return discretized delay distribution vector (if both input delays are
#'   vectors) or matrix
.convolve_two_delay_distributions <- function(first_delay,
                                          second_delay) {
  .are_valid_argument_values(list(
    list(first_delay, "delay_object", 1),
    list(second_delay, "delay_object", 1)
  )) #

  if (.is_numeric_vector(first_delay)) {
    if (.is_numeric_vector(second_delay)) {
      return(.convolve_delay_distribution_vectors(first_delay, second_delay))
    } else if (is.matrix(second_delay)) {
      return(.convolve_delay_distribution_vector_with_matrix(
        vector_a = first_delay,
        matrix_b = second_delay,
        vector_first = TRUE
      ))
    } else {
      stop("'second_delay' must be a numeric vector or matrix.")
    }
  } else if (is.matrix(first_delay)) {
    if (.is_numeric_vector(second_delay)) {
      return(.convolve_delay_distribution_vector_with_matrix(
        matrix_b = first_delay,
        vector_a = second_delay,
        vector_first = FALSE
      ))
    } else if (is.matrix(second_delay)) {
      # TODO work on .convolve_delay_distribution_matrices()
      stop("Convolution of two matrices is not available.")
      # return(.convolve_delay_distribution_matrices(first_delay, second_delay))
    } else {
      stop("'second_delay' must be a numeric vector or matrix.")
    }
  } else {
    stop("'first_delay' must be a numeric vector or matrix.")
  }
}

#' Convolve delay distributions
#'
#' Take a list of delay distributions and return their convolution.
#' The convolution of a delay A -> B and a delay B -> C corresponds to the
#' delay distribution of A -> C.
#' Delays are assumed to happen in the same chronological order as the order
#' they are given in in the \code{delays} list.
#'
#'
#' This function is flexible in the type of delay inputs it can handle.
#' Each delay in the \code{delays} list can be one of:
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
#' @example man/examples/convolve_delays.R
#' 
#' @inherit delay_high
#' @param n_report_time_steps integer. Length of incidence time series.
#' Use only when providing empirical delay data.
#' @inheritParams dating
#' @inheritDotParams get_matrix_from_empirical_delay_distr -empirical_delays -return_fitted_distribution
#' @inheritDotParams build_delay_distribution -distribution
#'
#' @return a discretized delay distribution vector or matrix.
#' A vector is returned when input delay distributions are constant through time:
#' either they are vectors already or in the form of a list-specified distribution.
#' A matrix is returned when at least one of the delays has a delay distribution
#' that can change through time. This is the case with empirical delay data
#' or if any of the input is already a delay distribution matrix.
#'
#' @export
convolve_delays <- function(delays,
                            n_report_time_steps = NULL,
                            ...) {

  .are_valid_argument_values(list(
    # We put '1' here, because we do not care here about checking the dimension of the matrix.
    list(delays, "delay_single_or_list", 1),
    list(n_report_time_steps, "null_or_int")
  ))


  dots_args <- .get_dots_as_list(...)

  if(.is_single_delay(delays)) {
    delay_distribution <- do.call(
      ".get_delay_distribution",
      c(
        list(
          delay = delays,
          n_report_time_steps = n_report_time_steps
        ),
        .get_shared_args(list(
          get_matrix_from_empirical_delay_distr,
          build_delay_distribution
        ), dots_args)
      )
    )
    return(delay_distribution)

  } else {

  delay_distributions <- lapply(delays,
                                function(delay) {
                                  delay_distribution <- do.call(
                                    ".get_delay_distribution",
                                    c(
                                      list(
                                        delay = delay,
                                        n_report_time_steps = n_report_time_steps
                                      ),
                                      .get_shared_args(list(
                                        get_matrix_from_empirical_delay_distr,
                                        build_delay_distribution
                                      ), dots_args)
                                    )
                                  )
                                })

  return(.convolve_delay_distributions(delay_distributions))
  }
}
