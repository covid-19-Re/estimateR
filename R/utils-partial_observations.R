#' Correct incidence data for yet-to-be-observed fraction of events
#'
#' Use this function to correct the tail of an incidence time series
#' if incidence was collected following a subsequent observation event.
#' For instance, if the incidence represents people starting to show symptoms of a disease
#' (dates of onset of symptoms), the data would typically have been collected among
#' individuals whose case was confirmed via a test.
#' If so, among all events of onset of symptoms, only those who had time to be
#' confirmed by a test were reported.
#' Thus, close to the present, there is an under-reporting of onset of symptoms events.
#' In order to account for this effect, this function divides each incidence value
#' by the probability of an event happening at a particular time step to have been observed.
#' Typically, this correction only affects the few most recent data points.
#' A trimming is done at the tail of the time series to avoid correcting for time steps
#' for which the observation probability is too low, which could result in too uncertain corrected values.
#' This trimming is tuned via the \code{cutoff_observation_probability} argument.
#'
#' The \code{ref_date} argument is only needed if the \code{delay_until_final_report}
#' is passed as a dataframe of individual delay observations (a.k.a empirical delay data).
#' In that case, \code{ref_date} must correspond to the date of the first time step in \code{incidence_data}.
#'
#' @inherit module_structure
#' @inherit delay_high
#'
#' @export
correct_for_partially_observed_data <- function(incidence_data,
                                                delay_until_final_report,
                                                cutoff_observation_probability = 0.1, # TODO set to a higher value i.e. 0.25 by default
                                                ref_date = NULL,
                                                time_step = "day",
                                                ...) {
  .are_valid_argument_values(list(
    list(incidence_data, "module_input"),
    list(delay_until_final_report, "delay_single_or_list", .get_input_length(incidence_data)),
    list(cutoff_observation_probability, "numeric_between_zero_one"),
    list(ref_date, "null_or_date"),
    list(time_step, "time_step")
  ))

  input <- .get_module_input(incidence_data)
  incidence_vector <- .get_values(input)

  dots_args <- .get_dots_as_list(...)

  delay_distribution_final_report <- do.call(
    "convolve_delays",
    c(
      list(
        delay = delay_until_final_report,
        n_report_time_steps = length(incidence_vector),
        ref_date = ref_date,
        time_step = time_step
      ),
      .get_shared_args(list(
        convolve_delays,
        build_delay_distribution,
        get_matrix_from_empirical_delay_distr), dots_args)
    )
  )

  if (NCOL(delay_distribution_final_report) == 1) {
    # delay_distribution_final_report is a vector, we build a delay distr matrix from it
    delay_distribution_matrix_final_report <- .get_matrix_from_single_delay_distr(delay_distribution_final_report,
      N = length(incidence_vector)
    )
  } else {
    # delay_distribution_final_report is a matrix, we truncate off the extra initial columns (required for R-L algo only)
    initial_offset <- ncol(delay_distribution_final_report) - length(incidence_vector) + 1
    delay_distribution_matrix_final_report <- delay_distribution_final_report[
      initial_offset:nrow(delay_distribution_final_report),
      initial_offset:ncol(delay_distribution_final_report)
    ]
  }

  Q_vector_observation_to_final_report <- apply(delay_distribution_matrix_final_report, MARGIN = 2, sum)

  # TODO improve this error
  # if (any(is.na(Q_vector_observation_to_final_report)) || isTRUE(any(Q_vector_observation_to_final_report == 0, na.rm = FALSE))) {
  if (any(is.na(Q_vector_observation_to_final_report))) {
    warning("Invalid delay_until_final_report argument.")
  }
  # TODO need to make sure that the matrix is the same size (as opposed to having extra columns leading)
  # TODO test with matrix delay
  # TODO we need to send an error if non zero incidence but zero probability of observation probably need a minimum value to replace zeroes by
  # TODO fix incidence vector if NaN or Inf values (there is cutoff anyway)
  incidence_vector <- incidence_vector / Q_vector_observation_to_final_report

  # Now we cut off values at the end of the time series,
  # those dates for which the probability of having observed an event that happened on that date is too low
  # We define 'too low' as being below a 'cutoff_observation_probability'
  tail_values_below_cutoff <- which(rev(Q_vector_observation_to_final_report) < cutoff_observation_probability)

  if (length(tail_values_below_cutoff) == 0) {
    cutoff <- 0
  } else {
    cutoff <- max(tail_values_below_cutoff)
  }

  truncated_incidence_vector <- incidence_vector[1:(length(incidence_vector) - cutoff)]

  return(.get_module_output(truncated_incidence_vector, input))
}
