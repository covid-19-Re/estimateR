#' Get relevant parameters from distribution
#'
#' @param f distribution function from stats package.
#' @inheritParams distribution
#'
#' @return list containing elements of \code{distribution} that overlap with arguments of \code{f}
.get_distribution_parms <- function(distribution, f) {
  # Remove the name element from the distribution list
  distribution <- within(distribution, rm("name"))

  # Only keep elements of 'distribution' that are arguments of function f
  distribution_parms <- distribution[names(distribution) %in% methods::formalArgs(f)]

  return(distribution_parms)
}

#' Get distribution function
#'
#' @param function_prefix character. 'd', 'q', 'p' or 'r'. see \code{\link[stats:Distributions]{stats::Distributions}}
#' @inheritParams distribution
#'
#' @return Density, distribution function, quantile function or random generation function for \code{distribution}
.get_distribution_function <- function(distribution, function_prefix) {
  .are_valid_argument_values(list(
    list(distribution, "distribution"),
    list(function_prefix, "function_prefix")
  ))
  f <- get(paste0(function_prefix, distribution[["name"]]), envir = loadNamespace("stats"))
  return(f)
}

#' Obtain quantile values for a distribution
#'
#' @inheritParams distribution
#' @param p vector of probabilities
#'
#' @return vector of quantiles
.get_quantiles <- function(distribution, p) {
  .are_valid_argument_values(list(
    list(distribution, "distribution"),
    list(p, "probability_distr_vector_high_tolerance")
  ))
  q_distribution_function <- .get_distribution_function(
    distribution = distribution,
    function_prefix = "q"
  )

  distribution_parms <- .get_distribution_parms(
    distribution = distribution,
    f = q_distribution_function
  )

  return(do.call(q_distribution_function, c(list(p = p), distribution_parms)))
}

#' Draw samples from a probability distribution.
#'
#' @inheritParams distribution
#' @param n integer. Number of samples to draw.
#'
#' @return vector containing \code{n} samples of \code{distribution}
.sample_from_distribution <- function(distribution, n) {
  .are_valid_argument_values(list(
    list(distribution, "distribution"),
    list(n, "integer")
  ))

  r_distribution_function <- .get_distribution_function(
    distribution = distribution,
    function_prefix = "r"
  )

  distribution_parms <- .get_distribution_parms(
    distribution = distribution,
    f = r_distribution_function
  )

  return(do.call(r_distribution_function, c(list(n = n), distribution_parms)))
}

#' Discretize a probability distribution.
#'
#' @param right_boundary positive numeric value.
#' Maximum number of time steps to discretize the \code{distribution} over.
#' @param offset_by_one boolean.
#' Set to TRUE if \code{distribution} represents the fit of data that was offset by one
#' (\code{fitted_data = original_data + 1}) to accommodate zeroes in \code{original_data}.
#' @inheritParams distribution
#'
#' @return vector containing weights of the discretized probability distribution.
.get_discretized_distribution <- function(distribution, right_boundary, offset_by_one) {
  .are_valid_argument_values(list(
    list(distribution, "distribution"),
    list(right_boundary, "positive_number"),
    list(offset_by_one, "boolean")
  ))

  p_distribution_function <- .get_distribution_function(
    distribution = distribution,
    function_prefix = "p"
  )

  distribution_parms <- .get_distribution_parms(
    f = p_distribution_function,
    distribution = distribution
  )

  if (offset_by_one) {
    x_values <- c(0, seq(from = 1.5, to = right_boundary, by = 1))
  } else {
    x_values <- c(0, seq(from = 0.5, to = right_boundary, by = 1))
  }

  cdf_values <- do.call(p_distribution_function, c(list(q = x_values), distribution_parms))

  if (length(cdf_values) == 1) {
    return(0)
  } else {
    return(diff(cdf_values))
  }
}

#' Get the number of steps before a quantile is reached.
#'
#' @inheritParams distribution
#' @param max_quantile numeric value. Upper quantile that needs to be reached.
#'
#' @return number of time steps required to for the probability distribution to reach \code{max_quantile}
.get_right_boundary_for_distribution_vector <- function(distribution, max_quantile) {
  .are_valid_argument_values(list(
    list(distribution, "distribution"),
    list(max_quantile, "numeric_between_zero_one")
  ))

  right_boundary <- ceiling(.get_quantiles(distribution, p = max_quantile)) + 1

  # Set the right boundary to at least two
  right_boundary <- max(right_boundary, 2)

  return(right_boundary)
}

#' Build a discretized probability distribution vector
#'
#' The discretization is done by integrating the probability density
#' on (0, 0.5), (0.5, 1.5), (1.5, 2.5)...
#'
#' @example man/examples/build_delay_distribution.R
#'
#' @inheritParams distribution
#' @inheritParams .get_discretized_distribution
#' @param max_quantile numeric value between 0 and 1.
#' Upper quantile reached by the last element in the discretized distribution vector.
#'
#' @return vector containing the discretized probability distribution vector of \code{distribution}.
#'
#' @export
build_delay_distribution <- function(distribution,
                                     max_quantile = 0.999,
                                     offset_by_one = FALSE) {
  .are_valid_argument_values(list(
    list(distribution, "distribution"),
    list(max_quantile, "numeric_between_zero_one"),
    list(offset_by_one, "boolean")
  ))

  right_boundary <- .get_right_boundary_for_distribution_vector(
    distribution = distribution,
    max_quantile = max_quantile
  )

  distribution_vector <- .get_discretized_distribution(
    distribution = distribution,
    right_boundary = right_boundary,
    offset_by_one = offset_by_one
  )

  return(distribution_vector)
}

#' Return probability distribution vector or matrix
#'
#' Can take a \code{distribution} list, a probability distribution vector,
#' a probability distribution matrix or empirical delay data as input.
#' If \code{delay} is already a delay distribution vector or matrix, it is returned as is.
#'
#' If \code{delay} is a \code{distribution} list,
#' this function builds and return the vector of discretized probability distribution.
#' If it is a vector, it checks that \code{delay}
#' is a valid discretized probability distribution and returns it.
#' Similarly, if \code{delay} is a matrix,
#' the function checks that it is in the correct format and returns it.
#' In a matrix, each column index corresponds to a time step associated with a date of event,
#' each row index corresponds to a time step associated with a date of event observation.
#' Each entry m_ij corresponds to the probability that an event occurring at time step j
#' is observed at time step i.
#' Matrices must be lower-triangular. Sums over columns must not exceed 1.
#'
#' See \code{\link{build_delay_distribution}} for details on the \code{distribution} list format;
#' see \code{\link{get_matrix_from_empirical_delay_distr}} for details on the empirical delay data format.
#'
#' @param delay list, vector, matrix or dataframe.
#' Delay distribution to transform or validate
#' into a vector of discretized probability distribution.
#' @inheritDotParams build_delay_distribution -distribution
#' @inherit dating
#' @inherit delay_empirical
#'
#' @return vector or matrix of discretized probability distribution.
.get_delay_distribution <- function(delay,
                                    n_report_time_steps = NULL,
                                    ref_date = NULL,
                                    time_step = "day",
                                    ...) {
  .are_valid_argument_values(list(
    list(delay, "delay_object", 1), # We put '1' here, because we do not care here about checking the dimension of the matrix.
    list(n_report_time_steps, "null_or_int"),
    list(ref_date, "null_or_date"),
    list(time_step, "time_step")
  ))

  dots_args <- .get_dots_as_list(...)

  if (is.data.frame(delay)) {
    if (is.null(n_report_time_steps) || n_report_time_steps == 0) {
      stop("Empirical delay data input but 'n_report_time_steps' parameter was not set or set to zero.")
    }

    delay_distribution <- do.call(
      "get_matrix_from_empirical_delay_distr",
      c(
        list(
          empirical_delays = delay,
          n_report_time_steps = n_report_time_steps,
          ref_date = ref_date,
          time_step = time_step
        ),
        .get_shared_args(list(get_matrix_from_empirical_delay_distr), dots_args)
      )
    )
  } else if (is.list(delay)) {
    delay_distribution <- do.call(
      "build_delay_distribution",
      c(
        list(distribution = delay),
        .get_shared_args(list(build_delay_distribution), dots_args)
      )
    )
  } else if (is.matrix(delay) || .is_numeric_vector(delay)) {
    delay_distribution <- delay
  } else {
    stop("Unknown delay type.")
  }

  return(delay_distribution)
}

#' Make delay distribution matrix from vector of delay distribution.
#'
#' @inheritParams distribution
#' @param N integer. Dimension of output matrix.
#'
#' @return discretized delay distribution matrix, representing a constant-through-time
#' delay distribution.
.get_matrix_from_single_delay_distr <- function(delay_distribution_vector, N) {
  .are_valid_argument_values(list(
    list(delay_distribution_vector, "probability_distr_vector"),
    list(N, "positive_integer")
  ))

  if (N >= length(delay_distribution_vector)) {
    delay_distribution_vector <- c(delay_distribution_vector, rep(0, times = N - length(delay_distribution_vector)))
  }

  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    delay_distribution_matrix[, i] <- c(rep(0, times = i - 1), delay_distribution_vector[1:(N - i + 1)])
  }

  return(delay_distribution_matrix)
}

# TODO maybe merge with .get_matrix_from_single_delay_distr by adding N parm and checking if list or unique vector
#' Build delay distribution matrix from list of delay distributions
#'
#' @param distributions list of distributions,
#' each element is either a distribution list or discretized probability distribution vector.
#' @inheritDotParams build_delay_distribution -distribution
#'
#' @return delay distribution matrix
.get_delay_matrix_from_delay_distribution_parms <- function(distributions, ...) {
  for (i in 1:length(distributions)) {
    .are_valid_argument_values(list(list(distributions[[i]], "distribution")))
  }
  dots_args <- .get_dots_as_list(...)

  # Generate list of delay distribution vectors
  delay_distribution_list <- lapply(distributions, function(distr) {
    do.call(
      "build_delay_distribution",
      c(
        list(distribution = distr),
        .get_shared_args(build_delay_distribution, dots_args)
      )
    )
  })

  N <- length(distributions)

  # Initialize empty matrix
  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)

  # Fill matrix by column
  for (i in 1:N) {
    delay_distr <- delay_distribution_list[[i]]

    # Right-pad delay_distr vector with zeroes if needed
    if (length(delay_distr) < N - i + 1) {
      delay_distr <- c(delay_distr, rep(0, times = N - i + 1 - length(delay_distr)))
    }
    delay_distribution_matrix[, i] <- c(rep(0, times = i - 1), delay_distr[1:(N - i + 1)])
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
                                             n_col_augment) {
  .are_valid_argument_values(list(
    list(delay_distribution_matrix, "probability_distr_matrix", 0),
    list(n_col_augment, "positive_integer")
  ))

  n_col_original <- ncol(delay_distribution_matrix)
  n_col_augmented <- n_col_original + n_col_augment

  # Initialize empty matrix
  augmented_matrix <- matrix(0, nrow = n_col_augmented, ncol = n_col_augmented)

  # Fill matrix by column

  # Start by duplicating first column in original matrix into 'n_col_augment' first columns of augmented_matrix
  for (i in 1:n_col_augment) {
    augmented_matrix[, i] <- c(rep(0, times = i - 1), delay_distribution_matrix[, 1], rep(0, times = n_col_augment - i + 1))
  }

  # Then fill with original matrix, adding the required zero on the top rows
  for (i in (n_col_augment + 1):n_col_augmented) {
    augmented_matrix[, i] <- c(rep(0, times = n_col_augment), delay_distribution_matrix[, i - n_col_augment])
  }

  return(augmented_matrix)
}

#' Get initial shift for deconvolution step
#'
#' Return the time step corresponding to a specified quantile
#' for a particular distribution (specified as a vector).
#'
#' @inheritParams distribution
#'
#' @return an integer value corresponding to the rounded value
#' corresponding to the specified quantile of the input delay distribution.
#' The returned calue is always non-negative.
.get_time_steps_quantile <- function(delay_distribution_vector, quantile = 0.5) {
  .are_valid_argument_values(
    list(
      list(delay_distribution_vector, "probability_distr_vector"),
      list(quantile, "numeric_between_zero_one")
    )
  )

  initial_shift <- ceiling(min(which(cumsum(delay_distribution_vector) > quantile))) - 1
  initial_shift <- max(initial_shift, 0, na.rm = TRUE)
  return(initial_shift)
}

#' Is input a single delay object?
#'
#' @param delay_list list or single delay object
#'
#' @return TRUE if single delay, FALSE otherwise
.is_single_delay <- function(delay_list) {
  .are_valid_argument_values(list(
    # We put '1' here, because we do not care here about checking the dimension of the matrix.
    list(delay_list, "delay_single_or_list", 1)
  ))

  if (is.list(delay_list) && !is.data.frame(delay_list)) {
    is_distribution <- try(.is_valid_distribution(delay_list, "dummy_name"), silent = TRUE)
    if ("try-error" %!in% class(is_distribution)) {
      return(TRUE)
    }
    return(FALSE)
  } else {
    return(TRUE)
  }
}
