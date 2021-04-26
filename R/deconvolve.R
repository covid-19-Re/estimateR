#TODO replace delay_incubation and delay_onset_to_report by list of delays

#TODO add parameter to allow delay_distribution to not sum up to 1 if not waiting time distribution (e.g. shedding load distribution)

#TODO redo doc
#' Infer Infection Events Dates from Delayed Observation
#'
#' This function reconstructs an incidence of infection events from incidence data representing delayed observations.
#' The assumption made is that delayed observations represent the convolution of the infections with a delay distribution.
#' \code{deconvolve_incidence} implements a deconvolution algorithm (Richardson-Lucy) to reconstruct
#' a vector of infection events from input data representing delayed observations.
#'
#'#TODO figure out input format
#'
#' @param incidence_data numeric.
#' @param deconvolution_method string. Options are "Richardson-Lucy delay distribution"
#' @param delay_incubation
#' @param delay_onset_to_report
#' @param simplify_output boolean. Return a numeric vector instead of module output object if output offset is zero.
#' @param ...
#' @param start_date
#' @param time_step
#' @param min_number_cases
#'
#' @return module output object.
#' @export
#'
#' @examples
#' #TODO add examples
deconvolve_incidence <- function( incidence_data,
                                  deconvolution_method = "Richardson-Lucy delay distribution",
                                  delay_incubation,
                                  delay_onset_to_report = c(1.0),
                                  start_date = NULL, 
                                  time_step = "day",
                                  min_number_cases = NULL,
                                  simplify_output = TRUE,
                                  ... ) {
  
  .are_valid_argument_values(list(list(user_input = incidence_data, input_type = "module_input", parameter_name = "incidence_data"),
                                  list(user_input=deconvolution_method, input_type="deconvolution_method", parameter_name="deconvolution_method"),
                                  list(user_input=delay_incubation, input_type="empirical_delay_data", parameter_name="delay_incubation"),
                                  list(user_input=delay_onset_to_report, input_type="empirical_delay_data", parameter_name="delay_onset_to_report"),
                                  list(user_input=start_date, input_type="null_or_date", parameter_name="start_date"),
                                  list(user_input=time_step, input_type="time_step", parameter_name="time_step"),
                                  list(user_input=min_number_cases, input_type="null_or_numeric", parameter_name="min_number_cases"),
                                  list(user_input=simplify_output, input_type="boolean", parameter_name="simplify_output")))
  
  
  input <- .get_module_input(incidence_data)

  #TODO generalize this to a list of delay inputs
  total_delay_distribution <- convolve_delay_inputs(delay_incubation,
                                                    delay_onset_to_report,
                                                    n_report_time_steps = .get_input_length(input),
                                                    start_date = start_date,
                                                    time_step = time_step,
                                                    min_number_cases = min_number_cases)

  if(deconvolution_method == "Richardson-Lucy delay distribution") {
    deconvolved_incidence <- .deconvolve_incidence_Richardson_Lucy(input,
                                                                   delay_distribution = total_delay_distribution,
                                                                   ... )
  } else {
    deconvolved_incidence <-  .make_empty_module_output()
  }

  if(simplify_output) {
    deconvolved_incidence <- .simplify_output(deconvolved_incidence)
  }

  return(deconvolved_incidence)
}

#' Deconvolve the incidence input with the Richardson-Lucy (R-L) algorithm
#'
#' @param incidence_input module input object.
#' @param delay_distribution numeric square matrix or vector.
#' @param threshold_chi_squared numeric. Threshold for chi-squared values under which the R-L algo stops.
#' @param max_iterations integer. Maximum threshold for the number of iterations in the Richardson-Lucy algo.
#' @param verbose Boolean. Print verbose output?
#'
#' @return module output object. Deconvolved incidence.
.deconvolve_incidence_Richardson_Lucy <- function(
  incidence_input,
  delay_distribution,
  threshold_chi_squared = 1,
  max_iterations = 100,
  verbose = FALSE
) {

  incidence_vector <- .get_values(incidence_input)

  # Test delay_distribution input
  #TODO apply proper validation step
  if(NCOL(delay_distribution) != 1 && NCOL(delay_distribution) != NROW(delay_distribution)) {
    #TODO cast proper error
    print("Delay distribution input must be a vector or square matrix")
  }

  length_original_vector <- length(incidence_vector)
  first_recorded_incidence <-  incidence_vector[1]
  last_recorded_incidence <- incidence_vector[length_original_vector]

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
