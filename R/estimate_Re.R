#TODO redoc with piecewise constant
#TODO add details about the two options for estimating Re
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
                         simplify_output = FALSE,
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
  } else if (estimation_method == "EpiEstim piecewise constant") {
    Re_estimate <- do.call(
      '.estimate_Re_EpiEstim_piecewise_constant',
      c(list(incidence_input = input),
        .get_shared_args(.estimate_Re_EpiEstim_piecewise_constant, dots_args))
    )
  } else {
    Re_estimate <- .make_empty_module_output()
  }

  if(simplify_output) {
    if(.is_list_of_outputs(Re_estimate)) {
    #TODO test if this works as intended
      Re_estimate <- merge_outputs(Re_estimate)
    } else {
      Re_estimate <- .simplify_output(Re_estimate)
    }
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
                        mean_Re_prior = 1,
                        output_HPD = FALSE) {

  .are_valid_argument_values(list(list(incidence_input, "module_input"),
                                  list(minimum_cumul_incidence, "non_negative_number"),
                                  list(estimation_window, "positive_integer"),
                                  list(mean_serial_interval, "positive_number"),
                                  list(std_serial_interval, "non_negative_number"),
                                  list(mean_Re_prior, "positive_number")))
  
  incidence_vector <- .get_values(incidence_input)

  if(sum(incidence_vector) < minimum_cumul_incidence) {
    stop("minimum_cumul_incidence parameter is set higher than total cumulative incidence.")
  }

  offset <- which(cumsum(incidence_vector) >= max(12, minimum_cumul_incidence))[1]
  offset <- max(2, estimation_window, ceiling(mean_serial_interval), offset)

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
  Re_estimate <- .get_module_output(R_instantaneous$R$`Mean(R)`,
                               incidence_input,
                              additional_offset)
  if(output_HPD){
    Re_highHPD <- .get_module_output(R_instantaneous$R$`Quantile.0.975(R)`,
                                 incidence_input,
                                 additional_offset)

    Re_lowHPD <- .get_module_output(R_instantaneous$R$`Quantile.0.025(R)`,
                                 incidence_input,
                                 additional_offset)

    return(list(Re_estimate=Re_estimate,
                Re_highHPD=Re_highHPD,
                Re_lowHPD=Re_lowHPD))
  } else {
    return(Re_estimate)
  }
}

#TODO doc
#' Estimate Re with EpiEstim in a piecewise-constant fashion
#'
#' @param minimum_cumul_incidence Numeric scalar. Minimum number of cumulated infections before starting the Re estimation
#' @param mean_serial_interval Numeric positive scalar. \code{mean_si} for \code{\link[EpiEstim]{estimate_R}}
#' @param std_serial_interval Numeric positive scalar. \code{std_si} for \code{\link[EpiEstim]{estimate_R}}
#' @param mean_Re_prior Numeric positive scalar. \code{mean prior} for \code{\link[EpiEstim]{estimate_R}}
#' @inheritParams inner_module
#'
#' @return module output object. mean of Re estimates
.estimate_Re_EpiEstim_piecewise_constant <- function(incidence_input,
                                                 minimum_cumul_incidence = 10,
                                                 interval_ends = NULL,
                                                 interval_length = 7,
                                                 mean_serial_interval = 4.8,
                                                 std_serial_interval  = 2.3,
                                                 mean_Re_prior = 1,
                                                 output_HPD = FALSE) {

  #TODO validate new args
  .are_valid_argument_values(list(list(incidence_input, "module_input"),
                                  list(minimum_cumul_incidence, "non_negative_number"),
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
  right_bound <- length(incidence_vector)

  if(!is.null(interval_ends)) {
    #TODO validate interval ends (not necessarily positive)
    # we make these be relative to index_offset (IMPORTANT, TODO doc)
    index_offset <- .get_offset(incidence_input)

    interval_ends <- sort(interval_ends - index_offset)

    interval_ends <- interval_ends[interval_ends > offset & interval_ends <= right_bound]
  } else {
    #TODO validate interval_length (must be positive integer)
    interval_ends <- seq(from = offset + interval_length -1, to = right_bound, by = interval_length)
    if(max(interval_ends) < right_bound) {
      interval_ends <- c(interval_ends, right_bound)
    }
  }

  if(length(interval_ends) < 1) {
    stop("No valid interval to estimate Re on.
         Check 'minimum_cumul_incidence', 'interval_ends' or 'interval_length' parameters.")
  }

  interval_starts <- c(offset, interval_ends[-length(interval_ends)] + 1)

  R_instantaneous <- EpiEstim::estimate_R(
    incidence_vector,
    method = "parametric_si",
    config = EpiEstim::make_config(
      list(
        mean_si = mean_serial_interval,
        std_si = std_serial_interval,
        t_start = interval_starts,
        t_end = interval_ends,
        mean_prior = mean_Re_prior)
    )
  )

  additional_offset <- interval_starts[1] - 1

  replicate_estimates_on_interval <- function(estimates, interval_starts, interval_ends) {
    replicated_estimates <- unlist(lapply(seq_along(interval_starts),
                                          function(x) {
                                            rep(estimates[x], interval_ends[x] - interval_starts[x] + 1)
                                          }))
    return(replicated_estimates)
  }

  Re_estimate <- replicate_estimates_on_interval(R_instantaneous$R$`Mean(R)`, interval_starts, interval_ends)
  Re_estimate <- .get_module_output(Re_estimate, incidence_input, additional_offset)

  if(output_HPD){

    Re_highHPD <- replicate_estimates_on_interval(R_instantaneous$R$`Quantile.0.975(R)`, interval_starts, interval_ends)
    Re_highHPD <- .get_module_output(Re_highHPD, incidence_input, additional_offset)

    Re_lowHPD <- replicate_estimates_on_interval(R_instantaneous$R$`Quantile.0.025(R)`, interval_starts, interval_ends)
    Re_lowHPD <- .get_module_output(Re_lowHPD, incidence_input, additional_offset)

    return(list(Re_estimate=Re_estimate,
                Re_highHPD=Re_highHPD,
                Re_lowHPD=Re_lowHPD))
  } else {
    return(Re_estimate)
  }
}
