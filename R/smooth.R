#TODO fill in documentation
#' Smooth Noisy Incidence Data
#'
#' @param incidence_data
#' @param smoothing_method
#'
#' @return
#' @export
#'
#' @examples
smooth_incidence <- function(incidence_data, smoothing_method = "LOESS", ...) {

  input <- get_module_input(incidence_data)

  #TODO have way to add additional parameters for each smoothing method
  smoothed_incidence <- dplyr::case_when(
    smoothing_method == "LOESS" ~ smooth_LOESS(input, ...),
    TRUE ~ rep(NA_real_, length.out = get_input_length(input))
  )
  return(smoothed_incidence)
}

#TODO fill in documentation
#' LOESS smoothing function
#'
#' @param incidence_input
#' @param days_incl integer
#' @param degree integer
#'
#' @return module output
#'
#' @examples
smooth_LOESS <- function(incidence_input, days_incl = 21, degree = 1) {
  incidence_vector <- incidence_input$values

  n_points <- length(incidence_vector)
  sel_span <- days_incl / n_points

  n_pad <- round(length(incidence_vector) * sel_span * 0.5)

  c_data <- data.frame(value = c(rep(0, n_pad), incidence_vector),
                       date_num = 1:(n_pad + n_points))

  c_data.lo <- stats::loess(value ~ date_num, data = c_data, span = sel_span, degree = degree)
  smoothed <- stats::predict(c_data.lo)
  smoothed[smoothed < 0] <- 0
  raw_smoothed_counts <- smoothed[(n_pad + 1):length(smoothed)]
  normalized_smoothed_counts <-
    raw_smoothed_counts * sum(incidence_vector, na.rm = T) / sum(raw_smoothed_counts, na.rm = T)

  return(get_module_output(normalized_smoothed_counts, incidence_input))
}
