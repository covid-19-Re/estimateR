#'Smooth noisy incidence data
#'
#'Currently only LOESS smoothing (\code{smoothing_method = "LOESS"}) is implemented.
#'
#'@inheritParams module_methods
#'@inheritParams module_structure
#'@inheritDotParams .smooth_LOESS -incidence_input
#'
#'@return module output. Smoothed incidence. TODO add details
#'@export
smooth_incidence <- function(incidence_data,
                             smoothing_method = "LOESS",
                             simplify_output = TRUE,
                             ...) {
  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(smoothing_method, "smoothing_method"),
                                  list(simplify_output, "boolean")))


  dots_args <- .get_dots_as_list(...)
  input <- .get_module_input(incidence_data)

  if(smoothing_method == "LOESS") {
    smoothed_incidence <- do.call(
      '.smooth_LOESS',
      c(list(incidence_input = input),
        .get_shared_args(.smooth_LOESS, dots_args))
    )
  } else {
    smoothed_incidence <- .make_empty_module_output()
  }

  if(simplify_output) {
    smoothed_incidence <- .simplify_output(smoothed_incidence)
  }

  return(smoothed_incidence)
}

#' LOESS smoothing function
#'
#' Prefer the use of the wrapper function \code{smooth_incidence(..., smoothing_method = "LOESS")}
#' instead of \code{.smooth_LOESS}.
#'
#' This function implements the LOESS method for smoothing noisy data.
#' It is essentially a wrapper around \code{\link[stats]{loess}}.
#' See the help section for \code{\link[stats]{loess}} for details on LOESS.
#'
#' #TODO add details on how data_points_incl relates to span.
#'
#' @inheritParams module_structure
#' @inheritParams inner_module
#' @param data_points_incl integer. Size of the window used in the LOESS algorithm.
#' @param degree integer. LOESS degree.
#'
#' @return module output. Smoothed incidence TODO add details.
.smooth_LOESS <- function(incidence_input, data_points_incl = 21, degree = 1) {

  .are_valid_argument_values(list(list(incidence_input, "module_input"),
                                  list(data_points_incl, "non_negative_number"), #minimal test; may have additional requirements
                                  list(degree, "non_negative_number"))) #minimal test; needs to be one of {0,1,2}, but stats::loess already throws if it isn't

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
