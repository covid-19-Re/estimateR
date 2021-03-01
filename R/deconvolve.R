
#TODO write documentation for each function
#TODO add a mandatory parameter for passing info on the delay distribution

#' Title
#'
#' @param incidence_data
#' @param estimation_method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
deconvolve_incidence <- function( incidence_data, deconvolution_method = "Richardson-Lucy delay distribution", ... ) {

  input <- get_module_input(incidence_data)

  deconvolved_incidence <- dplyr::case_when(
    deconvolution_method == "Richardson-Lucy delay distribution" ~ deconvolve_incidence_Richardson_Lucy(input, ... ),
    TRUE ~ rep(NA_real_, length.out = get_input_length(input))
  )
  return(deconvolved_incidence)
}

#TODO fill in doc
# incidence_vector is the vector with the original incidence
# delay_distribution_matrix is either a matrix or vector
#' Title
#'
#' @param incidence_input
#' @param delay_distribution
#' @param initial_delta
#' @param time_units_in_the_past
#' @param threshold_chi_squared
#' @param max_iterations
#' @param verbose
#'
#' @return
deconvolve_incidence_Richardson_Lucy <- function(
  incidence_input,
  delay_distribution,
  initial_delta,
  time_units_in_the_past = 30,
  threshold_chi_squared = 1,
  max_iterations = 100,
  verbose = FALSE
) {

  incidence_vector <- incidence_input$values


  #TODO test whether delay_distribution sums up to 1 if vector
  #TODO add parameter to allow delay_distribution to not sum up to 1 if not waiting time distribution

  # Test delay_distribution input
  if(NCOL(delay_distribution) != 1 && NCOL(delay_distribution) != NROW(delay_distribution)) {
    #TODO cast proper error
    print("Delay distribution input must be a vector or square matrix")
  }

  first_guess_delay <- ceiling(initial_delta)
  length_original_vector <- length(incidence_vector)
  first_recorded_incidence <-  incidence_vector[1]
  last_recorded_incidence <- incidence_vector[length(incidence_vector)]

  # prepare vector with initial guess for first step of deconvolution
  first_guess <- c(incidence_vector, rep(last_recorded_incidence, times = first_guess_delay))

  if(first_guess_delay < time_units_in_the_past) {
    first_guess <- c(rep(first_recorded_incidence, times = time_units_in_the_past - first_guess_delay), first_guess)
  } else if(time_units_in_the_past < first_guess_delay) {
    first_guess <- first_guess[-c(1:(first_guess_delay - time_units_in_the_past))]
  }

  original_incidence <- c(rep(0, times = length(first_guess) - length(incidence_vector)), incidence_vector)

  # Richardson-Lucy algorithm
  # initial step
  current_estimate <- first_guess
  chi_squared <- Inf
  count <- 1

  if(NCOL(delay_distribution) == 1) {
    delay_distribution_matrix <- get_matrix_from_constant_waiting_time_distr(delay_distribution,
                                                                             N=length(current_estimate))
  } else {
    delay_distribution_matrix <- delay_distribution[1:length(current_estimate), 1:length(current_estimate)]
  }

  truncated_delay_distribution_matrix <- delay_distribution_matrix[(1 + time_units_in_the_past):NROW(delay_distribution_matrix),, drop = F]

  Q_vector <- apply(truncated_delay_distribution_matrix, MARGIN = 2, sum)

  if (verbose) {
    cat("\tStart of Richardson-Lucy algorithm\n")
  }

  while(chi_squared > threshold_chi_squared & count <= max_iterations) {

    if (verbose) {
      cat("\t\tStep: ", count, " - Chi squared: ", chi_squared, "\n")
    }

    E <- as.vector(delay_distribution_matrix %*% current_estimate)
    B <- tidyr::replace_na(original_incidence/E, 0)

    current_estimate <- current_estimate / Q_vector *  as.vector(crossprod(B, delay_distribution_matrix))
    current_estimate <- tidyr::replace_na(current_estimate, 0)

    chi_squared <- 1/length_original_vector * sum((E[(time_units_in_the_past + 1): length(E)] - original_incidence[(time_units_in_the_past + 1) : length(original_incidence)])^2/E[(time_units_in_the_past + 1): length(E)], na.rm = T)
    count <- count + 1
  }

  additional_offset <- - min(first_guess_delay, time_units_in_the_past)
  # Remove last values as they cannot be properly inferred
  final_estimate <- current_estimate[1:(length(current_estimate) - first_guess_delay)]

  return(get_module_output(final_estimate, input, additional_offset))
}
