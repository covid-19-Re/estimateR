#' Estimate effective reproductive number Re from incidence data
#'
#' The incidence data input should represent infections,
#' as opposed to representing delayed observations of infections.
#' If the incidence data represents delayed observations of infections
#' then one should deconvolve it first with \code{deconvolve_incidence}.
#' TODO add details on inner function, that it wraps around EpiEstim and so on.
#'
#'
#' @inheritParams module_methods
#' @inheritParams module_structure
#' @inheritDotParams .estimate_Re_EpiEstim_sliding_window -incidence_input
#'
#' @return a module output object. Re estimates.
#' @export
estimate_Re <- function( incidence_data,
                         estimation_method = "EpiEstim sliding window",
                         simplify_output = TRUE,
                         ... ) {

  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(estimation_method, "estimation_method"),
                                  list(simplify_output, "boolean")))


  dots_args <- .get_dots_as_list(...)
  input <- .get_module_input(incidence_data)

  if(estimation_method == "EpiEstim sliding window") {
    Re_estimate <- do.call(
      '.estimate_Re_EpiEstim_sliding_window',
      c(list(incidence_input = input),
        .get_shared_args(.estimate_Re_EpiEstim_sliding_window, dots_args))
    )
  } else {
    Re_estimate <- .make_empty_module_output()
  }

  if(simplify_output) {
    Re_estimate <- .simplify_output(Re_estimate)
  }

  return(Re_estimate)
}


#TODO polish doc
#' Estimate Re with EpiEstim using a sliding window
#'
#' The Re value reported for time t corresponds to the value estimated
#' when assuming that is Re is constant over e.g. (T-3, T-2, T-1, T),
#' for a sliding window of 4 time steps.
#'
#' @param minimum_cumul_incidence Numeric scalar. Minimum number of cumulated infections before starting the Re estimation
#' @param estimation_window Integer scalar. Number of data points over which to assume Re to be constant.
#' @param mean_serial_interval Numeric positive scalar. \code{mean_si} for \code{\link[EpiEstim]{estimate_R}}
#' @param std_serial_interval Numeric positive scalar. \code{std_si} for \code{\link[EpiEstim]{estimate_R}}
#' @param mean_Re_prior Numeric positive scalar. \code{mean prior} for \code{\link[EpiEstim]{estimate_R}}
#' @inheritParams inner_module
#'
#' @return module output object. mean of Re estimates
.estimate_Re_EpiEstim_sliding_window <- function(incidence_input,
                        minimum_cumul_incidence = 10,
                        estimation_window = 3,
                        mean_serial_interval = 4.8,
                        std_serial_interval  = 2.3,
                        mean_Re_prior = 1) {

  .are_valid_argument_values(list(list(incidence_input, "module_input"),
                                  list(minimum_cumul_incidence, "non_negative_number"),
                                  list(estimation_window, "number"),
                                  list(mean_serial_interval, "number"),
                                  list(std_serial_interval, "non_negative_number"),
                                  list(mean_Re_prior, "number")))

  incidence_vector <- .get_values(incidence_input)

  if(sum(incidence_vector) < minimum_cumul_incidence) {
    stop("minimum_cumul_incidence parameter is set higher than total cumulative incidence.")
  }

  offset <- which(cumsum(incidence_vector) >= minimum_cumul_incidence)[1]
  # offset needs to be at least two for EpiEstim
  offset <- max(2, offset)

  right_bound <- length(incidence_vector) - (estimation_window - 1)

  if (is.na(offset) | right_bound < offset) {
    # No valid data point, return empty estimate
    return(c())
  }

  # Computation intervals corresponding to every position of the sliding window
  t_start <- seq(offset, right_bound)
  t_end <- t_start + estimation_window - 1

  R_instantaneous <- EpiEstim::estimate_R(
    incidence_vector,
    method = "parametric_si",
    config = EpiEstim::make_config(
      list(
        mean_si = mean_serial_interval,
        std_si = std_serial_interval,
        t_start = t_start,
        t_end = t_end,
        mean_prior = mean_Re_prior)
    )
  )

  additional_offset <- t_end[1] - 1
  R_mean <- .get_module_output(R_instantaneous$R$`Mean(R)`,
                               incidence_input,
                              additional_offset)

  return(R_mean)
}
