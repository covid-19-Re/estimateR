#' Infer infection events dates from delayed observations
#'
#' This function reconstructs an incidence of infection events
#' from incidence data representing delayed observations.
#' The assumption made is that delayed observations represent
#' the convolution of the time series of infections with a delay distribution.
#' \code{deconvolve_incidence} implements a deconvolution algorithm (Richardson-Lucy) to reconstruct
#' a vector of infection events from input data that represents delayed observations.
#'
#' @inheritParams module_methods
#' @inherit module_structure
#' @inheritParams delay_high
#' @inheritDotParams convolve_delays
#' @inheritDotParams .deconvolve_incidence_Richardson_Lucy -incidence_input
#'
#' @export
deconvolve_incidence <- function(incidence_data,
                                 deconvolution_method = "Richardson-Lucy delay distribution",
                                 delay,
                                 simplify_output = TRUE,
                                 ...) {

  # TODO be careful if we relax the constraints on incidence_data:
  # .get_input_length(incidence_data) may not make sense anymore
  # TODO need to make sure whether delays can be matrices with nrows > length(incidence_data)
  .are_valid_argument_values(list(
    list(incidence_data, "module_input"),
    list(deconvolution_method, "deconvolution_method"),
    list(simplify_output, "boolean"),
    list(delay, "delay_single_or_list", .get_input_length(incidence_data))
  ))

  dots_args <- .get_dots_as_list(...)
  input <- .get_module_input(incidence_data)

  total_delay_distribution <- do.call(
    "convolve_delays",
    c(
      list(
        delays = delay,
        n_report_time_steps = .get_input_length(input)
      ),
      .get_shared_args(list(
        convolve_delays,
        build_delay_distribution,
        get_matrix_from_empirical_delay_distr), dots_args)
    )
  )

  if (deconvolution_method == "Richardson-Lucy delay distribution") {
    deconvolved_incidence <- do.call(
      ".deconvolve_incidence_Richardson_Lucy",
      c(
        list(
          incidence_input = input,
          delay_distribution = total_delay_distribution
        ),
        .get_shared_args(.deconvolve_incidence_Richardson_Lucy, dots_args)
      )
    )
  } else {
    deconvolved_incidence <- .make_empty_module_output()
  }

  if (simplify_output) {
    deconvolved_incidence <- .simplify_output(deconvolved_incidence)
  }

  return(deconvolved_incidence)
}

#' Deconvolve the incidence input with the Richardson-Lucy (R-L) algorithm
#'
#' @inheritParams inner_module
#' @inheritParams universal_params
#' @param delay_distribution numeric square matrix or vector.
#' @param threshold_chi_squared numeric scalar. Threshold for chi-squared values under which the R-L algorithm stops.
#' @param max_iterations integer. Maximum threshold for the number of iterations in the R-L algorithm.
#' @inherit module_structure
.deconvolve_incidence_Richardson_Lucy <- function(incidence_input,
                                                  delay_distribution,
                                                  threshold_chi_squared = 1,
                                                  max_iterations = 100,
                                                  verbose = FALSE) {

  .are_valid_argument_values(list(
    list(incidence_input, "module_input"),
    list(delay_distribution, "computation_ready_delay_object", .get_input_length(incidence_input)),
    list(threshold_chi_squared, "non_negative_number"),
    list(max_iterations, "non_negative_number"),
    list(verbose, "boolean")
  ))

  incidence_vector <- .get_values(incidence_input)
  length_original_vector <- length(incidence_vector)
  first_recorded_incidence <- incidence_vector[1]
  last_recorded_incidence <- incidence_vector[length_original_vector]

  if (NCOL(delay_distribution) == 1) { # delay_distribution is not a matrix yet.
    n_time_units_left_extension <- .get_time_steps_quantile(delay_distribution, quantile = 0.99)
    initial_shift <- .get_time_steps_quantile(delay_distribution, quantile = 0.5)

    delay_distribution_matrix <- .get_matrix_from_single_delay_distr(
      delay_distribution,
      N = length(incidence_vector) + n_time_units_left_extension
    )
  } else {
    delay_distribution_matrix <- delay_distribution

    if (NCOL(delay_distribution_matrix) < length_original_vector) {
      stop("The dimension of 'delay_distribution' cannot be smaller than the length of 'incidence_input'.")
    }
    n_time_units_left_extension <- NCOL(delay_distribution_matrix) - length_original_vector

    initial_shift <- min(
      n_time_units_left_extension,
      .get_time_steps_quantile(delay_distribution_matrix[, 1], quantile = 0.5)
    )
  }

  # Here we could decide to either extend with zeroes when we know it's zero, or with an extrapolation of the early values
  # With zeroes does it have the benefit that whatever value we have in the deconvolved values, there is no effect on optim step?
  original_incidence <- c(rep(0, times = n_time_units_left_extension), incidence_vector)

  ### Richardson-Lucy algorithm

  ## Initial step
  # Prepare vector with initial guess for first step of deconvolution
  # Here we could also decide to extend with extrapolation of last values
  current_estimate <- c(incidence_vector, rep(last_recorded_incidence, times = initial_shift))

  extra_left_steps <- n_time_units_left_extension - initial_shift

  if (extra_left_steps > 0) {
    current_estimate <- c(rep(first_recorded_incidence, times = extra_left_steps), current_estimate)
  } else if (extra_left_steps < 0) {
    stop("Initial shift in R-L algo should be less than the number of steps padded on the left side.")
  }

  chi_squared <- Inf
  count <- 1

  truncated_delay_distribution_matrix <- delay_distribution_matrix[(1 + n_time_units_left_extension):NROW(delay_distribution_matrix), , drop = F]

  observation_probability <- apply(truncated_delay_distribution_matrix, MARGIN = 2, sum)

  if (verbose) {
    cat("\tStart of Richardson-Lucy algorithm\n")
  }

  ## Iterative steps
  while (chi_squared > threshold_chi_squared & count <= max_iterations) {
    if (verbose) {
      cat("\t\tStep: ", count, " - Chi squared: ", chi_squared, "\n")
    }

    convolved_estimate <- as.vector(delay_distribution_matrix %*% current_estimate)
    ratio_original_to_reconstructed <- tidyr::replace_na(original_incidence / convolved_estimate, 0)

    current_estimate <- current_estimate / observation_probability * as.vector(crossprod(ratio_original_to_reconstructed, delay_distribution_matrix))
    current_estimate <- tidyr::replace_na(current_estimate, 0)

    chi_squared <- 1 / length_original_vector *
      sum((convolved_estimate[(n_time_units_left_extension + 1):length(convolved_estimate)] -
             original_incidence[(n_time_units_left_extension + 1):length(original_incidence)])^2 /
            convolved_estimate[(n_time_units_left_extension + 1):length(convolved_estimate)], na.rm = T)

    count <- count + 1
  }

  if (verbose) {
    cat("\tEnd of Richardson-Lucy algorithm\n")
  }

  additional_offset <- -initial_shift
  # Remove first and last values as they cannot be properly inferred
  final_estimate <- current_estimate[(1 + extra_left_steps):(length(current_estimate) - initial_shift)]

  return(.get_module_output(final_estimate,.get_offset(incidence_input), additional_offset))
}
