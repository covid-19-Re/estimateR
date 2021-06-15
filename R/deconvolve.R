#' Infer Infection Events Dates from Delayed Observation
#'
#' This function reconstructs an incidence of infection events from incidence data representing delayed observations.
#' The assumption made is that delayed observations represent the convolution of the infections with a delay distribution.
#' \code{deconvolve_incidence} implements a deconvolution algorithm (Richardson-Lucy) to reconstruct
#' a vector of infection events from input data representing delayed observations.
#'
#' @inheritParams module_methods
#' @inheritParams module_structure
#' @inheritParams delay_high
#' @inheritDotParams convolve_delay_inputs
#' @inheritDotParams .deconvolve_incidence_Richardson_Lucy -incidence_input
#'
#' @return module output object. Inferred incidence of infection events.
#' @export
deconvolve_incidence <- function( incidence_data,
                                  deconvolution_method = "Richardson-Lucy delay distribution",
                                  delay_incubation,
                                  delay_onset_to_report = c(1.0),
                                  simplify_output = TRUE,
                                  ... ) {

  #TODO be careful if we relax the constraints on incidence_data:
  #.get_input_length(incidence_data) make not make sense anymore
  #TODO need to make sure whether delays can be matrices with nrows > length(incidence_data)
  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(deconvolution_method, "deconvolution_method"),
                                  list(delay_incubation, "delay_object", .get_input_length(incidence_data)),
                                  list(delay_onset_to_report, "delay_object", .get_input_length(incidence_data)),
                                  list(simplify_output, "boolean")))



  dots_args <- .get_dots_as_list(...)
  input <- .get_module_input(incidence_data)

  #TODO skip convolution if only one delay
  total_delay_distribution <- do.call(
    'convolve_delay_inputs',
    c(list(delay_incubation = delay_incubation,
           delay_onset_to_report = delay_onset_to_report,
           n_report_time_steps = .get_input_length(input)),
      .get_shared_args(convolve_delay_inputs, dots_args))
  )

  if(deconvolution_method == "Richardson-Lucy delay distribution") {

    deconvolved_incidence <- do.call(
      '.deconvolve_incidence_Richardson_Lucy',
      c(list(incidence_input = input,
             delay_distribution = total_delay_distribution),
        .get_shared_args(.deconvolve_incidence_Richardson_Lucy, dots_args))
    )
  } else {
    deconvolved_incidence <-  .make_empty_module_output()
  }

  if(simplify_output) {
    deconvolved_incidence <- .simplify_output(deconvolved_incidence)
  }

  return(deconvolved_incidence)
}

#TODO test estimates with onset dates.
#TODO rework on verbosity
#' Deconvolve the incidence input with the Richardson-Lucy (R-L) algorithm
#'
#' @inheritParams inner_module
#' @inheritParams universal_params
#' @param delay_distribution numeric square matrix or vector. TODO refactor to estimateR
#' @param threshold_chi_squared numeric scalar. Threshold for chi-squared values under which the R-L algorithm stops.
#' @param max_iterations integer. Maximum threshold for the number of iterations in the R-L algorithm.
#'
#' @return module output object. Deconvolved incidence.
.deconvolve_incidence_Richardson_Lucy <- function(
  incidence_input,
  delay_distribution,
  threshold_chi_squared = 1,
  max_iterations = 100,
  verbose = FALSE
) {

  #TODO delay_distribution must be vector or matrix (not delay data or distribution list)
  .are_valid_argument_values(list(list(incidence_input, "module_input"),
                                  list(delay_distribution, "delay_object", .get_input_length(incidence_input)),
                                  list(threshold_chi_squared, "non_negative_number"),
                                  list(max_iterations, "non_negative_number"),
                                  list(verbose, "boolean")))

  incidence_vector <- .get_values(incidence_input)

  # Test delay_distribution input
  #TODO apply proper validation step
  if(NCOL(delay_distribution) != 1 && NCOL(delay_distribution) != NROW(delay_distribution)) {
    stop("Delay distribution input must be a vector or square matrix")
  }

  length_original_vector <- length(incidence_vector)
  first_recorded_incidence <-  incidence_vector[1]
  last_recorded_incidence <- incidence_vector[length_original_vector]

  #TODO reorganize this: build matrix beforehand and get shift uniquely from matrix size
  if(NCOL(delay_distribution) == 1) {
    first_guess_delay <- .get_initial_deconvolution_shift(delay_distribution)
  } else {
    #TODO add check that NCOL(delay_distribution) > length_original_vector
    first_guess_delay <- NCOL(delay_distribution) - length_original_vector
  }

  original_incidence <- c(rep(0, times = first_guess_delay), incidence_vector)

  ### Richardson-Lucy algorithm
  #TODO document the math notations used in the algo and possibly rename some intermediary variables (E, B...)

  ## Initial step
  # Prepare vector with initial guess for first step of deconvolution
  current_estimate <- c(incidence_vector, rep(last_recorded_incidence, times = first_guess_delay))
  chi_squared <- Inf
  count <- 1

  if(NCOL(delay_distribution) == 1) {
    delay_distribution_matrix <- .get_matrix_from_single_delay_distr(delay_distribution,
                                                                     N=length(current_estimate))
  } else {
    delay_distribution_matrix <- delay_distribution
  }

  truncated_delay_distribution_matrix <- delay_distribution_matrix[(1 + first_guess_delay):NROW(delay_distribution_matrix),, drop = F]

  Q_vector <- apply(truncated_delay_distribution_matrix, MARGIN = 2, sum)

  if (verbose) {
    cat("\tStart of Richardson-Lucy algorithm\n")
  }

  ## Iterative steps
  while(chi_squared > threshold_chi_squared & count <= max_iterations) {

    if (verbose) {
      cat("\t\tStep: ", count, " - Chi squared: ", chi_squared, "\n")
    }

    E <- as.vector(delay_distribution_matrix %*% current_estimate)
    B <- tidyr::replace_na(original_incidence/E, 0)

    current_estimate <- current_estimate / Q_vector *  as.vector(crossprod(B, delay_distribution_matrix))
    current_estimate <- tidyr::replace_na(current_estimate, 0)

    chi_squared <- 1/length_original_vector * sum((E[(first_guess_delay + 1): length(E)] - original_incidence[(first_guess_delay + 1) : length(original_incidence)])^2/E[(first_guess_delay + 1): length(E)], na.rm = T)
    count <- count + 1
  }

  if (verbose) {
    cat("\tEnd of Richardson-Lucy algorithm\n")
  }

  additional_offset <- - first_guess_delay
  # Remove last values as they cannot be properly inferred
  final_estimate <- current_estimate[1:(length(current_estimate) - first_guess_delay)]

  return(.get_module_output(final_estimate, incidence_input, additional_offset))
}
