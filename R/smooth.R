#TODO figure out how to deal with sometimes having incidence that needs to be left-padded with zeroes and sometimes not

#' Smooth Noisy Incidence Data
#'
#' Currently only LOESS smoothing (smoothing_method = "LOESS")
#' is implemented.
#'
#'#TODO specify input format
#' @param incidence_data
#' @param smoothing_method string. Options are "LOESS".
#'
#' @return module output. Smoothed incidence.
#' @export
#'
#'#TODO fill in example
#' @examples
smooth_incidence <- function(incidence_data, smoothing_method = "LOESS", ...) {

  input <- .get_module_input(incidence_data)

  smoothed_incidence <- dplyr::case_when(
    smoothing_method == "LOESS" ~ .smooth_LOESS(input, ...),
    TRUE ~ rep(NA_real_, length.out = .get_input_length(input))
  )
  return(smoothed_incidence)
}

#' LOESS smoothing function
#'
#' see help page of \link{stats::loess} for more details.
#'
#' #TODO add details on how data_points_incl relates to span.
#'
#' @param incidence_input module input. Noisy incidence data to smooth.
#' @param data_points_incl integer. Size of the window used in the LOESS algorithm.
#' @param degree integer. LOESS degree.
#'
#' @return module output. Smoothed incidence.
.smooth_LOESS <- function(incidence_input, data_points_incl = 21, degree = 1) {
  incidence_vector <- .get_values(incidence_input)

  n_points <- length(incidence_vector)
  sel_span <- data_points_incl / n_points

  n_pad <- round(length(incidence_vector) * sel_span * 0.5)

  c_data <- data.frame(value = c(rep(0, n_pad), incidence_vector),
                       date_num = 1:(n_pad + n_points))

  c_data.lo <- stats::loess(value ~ date_num, data = c_data, span = sel_span, degree = degree)
  smoothed <- stats::predict(c_data.lo)
  smoothed[smoothed < 0] <- 0
  raw_smoothed_counts <- smoothed[(n_pad + 1):length(smoothed)]
  normalized_smoothed_counts <-
    raw_smoothed_counts * sum(incidence_vector, na.rm = T) / sum(raw_smoothed_counts, na.rm = T)

  return(.get_module_output(normalized_smoothed_counts, incidence_input))
}
