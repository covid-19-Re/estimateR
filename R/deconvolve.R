#' Infer Infection Events Dates from Delayed Observation
#'
#' This function reconstructs an incidence of infection events from incidence data representing delayed observations.
#' The assumption made is that delayed observations represent the convolution of the infections with a delay distribution.
#' \code{deconvolve_incidence} implements a deconvolution algorithm (Richardson-Lucy) to reconstruct
#' a vector of infection events from input data representing delayed observations.
#'
#'#TODO figure out input format
#' @param incidence_data numeric.
#' @param deconvolution_method string. Options are "Richardson-Lucy delay distribution"
#' @param simplify_output boolean. Return a numeric vector instead of module output object if output offset is zero.
#' @param ...
#'
#' @return module output object.
#' @export
#'
#' @examples
#' #TODO add examples
#TODO add a mandatory parameter for passing info on the delay distribution
deconvolve_incidence <- function( incidence_data, deconvolution_method = "Richardson-Lucy delay distribution", simplify_output = TRUE, ... ) {

  input <- .get_module_input(incidence_data)

  if(deconvolution_method == "Richardson-Lucy delay distribution") {
    deconvolved_incidence <- .deconvolve_incidence_Richardson_Lucy(input, ... )
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

  #TODO test whether delay_distribution sums up to 1 if vector
  #TODO add parameter to allow delay_distribution to not sum up to 1 if not waiting time distribution
  #TODO compute first_guess_delay if delay_distribution is a matrix
  ##create a general utility function that does it for both vectors and matrices, with customisable way of getting the delay (mode, median, mean,..)

  # Test delay_distribution input
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
