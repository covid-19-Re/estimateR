#TODO replace mean/std serial interval inputs by a single input combining the two

#' Estimate effective reproductive number Re from incidence data
#'
#' The incidence data input should represent infections,
#' as opposed to representing delayed observations of infections.
#' If the incidence data represents delayed observations of infections
#' then one should deconvolve it first with \code{deconvolve_incidence}.
#'
#' For now, only one Re estimation function is implemented.
#'
#' #TODO specify input format
#'
#' @param incidence_data Format still to define.
#' @param estimation_method string. Options are "EpiEstim sliding window".
#' @param simplify_output boolean. Return a numeric vector instead of module output object if output offset is zero.
#' @param ...
#'
#'#TODO figure out output format
#' @return a module output object (subject to change)
#' @export
#'
#' @examples
#' #TODO add examples
estimate_Re <- function( incidence_data, estimation_method = "EpiEstim sliding window", simplify_output = TRUE, ... ) {

  .are_valid_argument_values(as.list(environment()))
  input <- .get_module_input(incidence_data)

  if(estimation_method == "EpiEstim sliding window") {
    Re_estimate <- .estimate_Re_EpiEstim_sliding_window(input, ... )
  } else {
    Re_estimate <- .make_empty_module_output()
  }

  if(simplify_output) {
    Re_estimate <- .simplify_output(Re_estimate)
  }

  return(Re_estimate)
}


#' Estimate Re with EpiEstim using a sliding window
#'
#' The Re value reported for time t corresponds to the value estimated
#' when assuming that is Re is constant over e.g. (T-3, T-2, T-1, T),
#' for a sliding window of 4 time steps.
#'
#' @param incidence_input module input object. Time series of infections.
#' @param minimul_cumul_incidence Numeric scalar. Minimum number of cumulated infections before starting the Re estimation
#' @param estimation_window Integer scalar. Number of data points over which to assume Re to be constant.
#' @param mean_serial_interval Numeric positive scalar. "mean_si" for EpiEstim::estimate_R
#' @param std_serial_interval Numeric positive scalar. "std_si" for EpiEstim::estimate_R
#' @param mean_Re_prior Numeric positive scalar. "mean prior" for EpiEstim::estimate_R
#'
#' @return module output object. mean of Re estimates
.estimate_Re_EpiEstim_sliding_window <- function(incidence_input,
                        minimul_cumul_incidence = 0,
                        estimation_window = 3,
                        mean_serial_interval = 4.8,
                        std_serial_interval  = 2.3,
                        mean_Re_prior = 1) {

  incidence_vector <- .get_values(incidence_input)

  # Ensure that incidence_vector has no NA or negative values
  #TODO make these checks in the get_module_input function
  incidence_vector <- dplyr::if_else(is.na(incidence_vector) | incidence_vector < 0, 0, incidence_vector)

  offset <- which(cumsum(incidence_vector) >= minimul_cumul_incidence)[1]
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
