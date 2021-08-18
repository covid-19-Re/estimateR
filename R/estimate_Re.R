#' Estimate the effective reproductive number Re through time from incidence data
#'
#' \code{estimate_Re()} takes the number of infections through time
#' and computes the Re value through time (also known as Rt).
#'
#' The incidence input should represent infections,
#' as opposed to representing delayed observations of infections.
#' If the incidence data represents delayed observations of infections,
#' one should first reconstruct the incidence of infections using
#' \code{deconvolve_incidence()} or \code{get_infections_from_incidence()}
#' which wraps around it and includes a smoothing step of the delayed observations.
#'
#' \code{estimate_Re()} wraps around the \code{estimate_R()} function
#' of the \code{EpiEstim} package from Cori et al, 2013.
#' \code{estimate_Re()} allows for two types of Re estimations:
#' \enumerate{
#' \item A sliding-window estimation.
#' For each time step T, the Re(T) value is computed by assuming that Re is constant
#' for time steps (T-X+1, ..., T-1, T), with X being the sliding window size.
#' This option is chosen by setting \code{estimation_method = "EpiEstim sliding window"}
#' \item A piecewise-constant estimation.
#' Re(t) is computed as being a piecewise-constant function of time.
#' The length of each step can be a fixed number of time steps.
#' That number is specified using the \code{interval_length} parameter.
#' The length of each step can also be irregular.
#' This can be useful if the boundaries of the steps are meant to coincide with particular events
#' such as the implementation or cancellation of public health interventions.
#' The right boundaries of the steps are specified using the \code{interval_ends} parameter.
#' This option is chosen by setting \code{estimation_method = "EpiEstim piecewise constant"}
#' }
#'
#'
#' @param simplify_output boolean. Simplify the output when possible?
#' @inheritParams module_methods
#' @inheritParams module_structure
#' @inheritDotParams .estimate_Re_EpiEstim_sliding_window -incidence_input
#' @inheritDotParams .estimate_Re_EpiEstim_piecewise_constant -incidence_input
#'
#' @return A list with two elements:
#'  \enumerate{
#'  \item A numeric vector named \code{values}: the result of the computations on the input data.
#'  \item An integer named \code{index_offset}: the offset, counted in number of time steps,
#'  by which the result is shifted compared to an \code{index_offset} of \code{0}.
#'  This parameter allows one to keep track of the date of the first value in \code{values}
#'  without needing to carry a \code{date} column around.
#'  A positive offset means \code{values} are delayed in the future compared to the reference values.
#'  A negative offset means the opposite.
#'  Note that the \code{index_offset} of the output of the function call
#'  accounts for the (optional) \code{index_offset} of the input.
#'  }
#'  If \code{index_offset} is \code{0} and \code{simplify_output = TRUE},
#'  the \code{index_offset} is dropped and the \code{values}
#'  element is returned as a numeric vector.
#'
#'  If \code{output_HPD = TRUE} (additional parameter),
#'  the highest posterior density interval boundaries are output along with the mean Re estimates.
#'  In that case, a list of three lists is returned:
#' \itemize{
#' \item \code{Re_estimate} contains the Re estimates.
#' \item \code{Re_highHPD} and \code{Re_lowHPD} contain
#' the higher and lower boundaries of the HPD interval,
#' as computed by \code{\link[EpiEstim]{estimate_R}}
#' }
#' If, in addition, \code{simplify_output = TRUE},
#' then the 3 elements are merged into a single dataframe by \code{merge_outputs()}.
#' A date column can be added to the dataframe by passing an extra \code{ref_date} argument
#' (see \code{\link{merge_outputs}} for details).
#'
#'
#' @export
estimate_Re <- function(incidence_data,
                        estimation_method = "EpiEstim sliding window",
                        simplify_output = FALSE,
                        ...) {
  .are_valid_argument_values(list(
    list(incidence_data, "module_input"),
    list(estimation_method, "estimation_method"),
    list(simplify_output, "boolean")
  ))


  dots_args <- .get_dots_as_list(...)
  input <- .get_module_input(incidence_data)

  if (estimation_method == "EpiEstim sliding window") {
    Re_estimate <- do.call(
      ".estimate_Re_EpiEstim_sliding_window",
      c(
        list(incidence_input = input),
        .get_shared_args(.estimate_Re_EpiEstim_sliding_window, dots_args)
      )
    )
  } else if (estimation_method == "EpiEstim piecewise constant") {
    Re_estimate <- do.call(
      ".estimate_Re_EpiEstim_piecewise_constant",
      c(
        list(incidence_input = input),
        .get_shared_args(.estimate_Re_EpiEstim_piecewise_constant, dots_args)
      )
    )
  } else {
    Re_estimate <- .make_empty_module_output()
  }

  if (simplify_output) {
    if (.is_list_of_outputs(Re_estimate)) {
      Re_estimate <- do.call(
        "merge_outputs",
        c(list(output_list = Re_estimate),
          .get_shared_args(merge_outputs, dots_args)
        )
      )
    } else {
      Re_estimate <- .simplify_output(Re_estimate)
    }
  }

  return(Re_estimate)
}

#' Estimate Re with EpiEstim using a sliding window
#'
#' The Re value reported for time t corresponds to the value estimated
#' when assuming that is Re is constant over e.g. (T-3, T-2, T-1, T),
#' for a sliding window of 4 time steps.
#'
#' @param estimation_window Use with \code{estimation_method = "EpiEstim sliding window"}
#' Positive integer value.
#' Number of data points over which to assume Re to be constant.
#' @inheritParams inner_module
#' @inherit EpiEstim_wrapper
.estimate_Re_EpiEstim_sliding_window <- function(incidence_input,
                                                 import_incidence_input = NULL,
                                                 minimum_cumul_incidence = 12,
                                                 estimation_window = 3,
                                                 mean_serial_interval = 4.8,
                                                 std_serial_interval = 2.3,
                                                 mean_Re_prior = 1,
                                                 output_HPD = FALSE) {
  .are_valid_argument_values(list(
    list(incidence_input, "module_input"),
    list(minimum_cumul_incidence, "non_negative_number"),
    list(estimation_window, "positive_integer"),
    list(mean_serial_interval, "positive_number"),
    list(std_serial_interval, "non_negative_number"),
    list(mean_Re_prior, "positive_number")
  ))

  incidence_vector <- .get_values(incidence_input)
  if (sum(incidence_vector) < minimum_cumul_incidence) {
    stop("minimum_cumul_incidence parameter is set higher than total cumulative incidence.")
  }

  if(!is.null(import_incidence_input)) {
    .are_valid_argument_values(list(
      list(import_incidence_input, "module_input")
    ))

    incidence <- merge_outputs(output_list = list(local = incidence_input,
                                imported = import_incidence_input),
                           include_index = FALSE) %>%
      tidyr::replace_na(list(local = 0, imported = 0))

    incidence_length <- nrow(incidence)
    input_offset <- min(.get_offset(incidence_input), .get_offset(import_incidence_input))

  } else {
    incidence <- incidence_vector
    incidence_length <- length(incidence)
    input_offset <- .get_offset(incidence_input)
  }

  offset <- which(cumsum(incidence_vector) >= minimum_cumul_incidence)[1]
  # We use the criteria on the offset from Cori et al. 2013
  # (and offset needs to be at least two for EpiEstim)
  offset <- max(estimation_window, ceiling(mean_serial_interval), offset, 2)

  right_bound <- incidence_length - (estimation_window - 1)

  if (is.na(offset) | right_bound < offset) {
    # No valid data point, return empty estimate
    return(c())
  }

  # Computation intervals corresponding to every position of the sliding window
  t_start <- seq(offset, right_bound)
  t_end <- t_start + estimation_window - 1

  R_instantaneous <- suppressWarnings(EpiEstim::estimate_R(
    incid = incidence,
    method = "parametric_si",
    config = EpiEstim::make_config(
      list(
        mean_si = mean_serial_interval,
        std_si = std_serial_interval,
        t_start = t_start,
        t_end = t_end,
        mean_prior = mean_Re_prior
      )
    )
  )
  )

  additional_offset <- t_end[1] - 1
  Re_estimate <- .get_module_output(
    R_instantaneous$R$`Mean(R)`,
    input_offset,
    additional_offset
  )
  if (output_HPD) {
    Re_highHPD <- .get_module_output(
      R_instantaneous$R$`Quantile.0.975(R)`,
      input_offset,
      additional_offset
    )

    Re_lowHPD <- .get_module_output(
      R_instantaneous$R$`Quantile.0.025(R)`,
      input_offset,
      additional_offset
    )

    return(list(
      Re_estimate = Re_estimate,
      Re_highHPD = Re_highHPD,
      Re_lowHPD = Re_lowHPD
    ))
  } else {
    return(Re_estimate)
  }
}

#' Estimate Re with EpiEstim in a piecewise-constant fashion
#'
#' This function returns piecewise-constant Re estimates.
#'
#' @param interval_ends Use with \code{estimation_method = "EpiEstim piecewise constant"}
#' Integer vector. Optional argument.
#' If provided, \code{interval_ends} overrides the \code{interval_length} argument.
#' Each element of \code{interval_ends} specifies the right boundary
#' of an interval over which Re is assumed to be constant for the calculation.
#' Values in \code{interval_ends} must be integer values corresponding
#' with the same numbering of time steps as given by \code{incidence_input}.
#' In other words, \code{interval_ends} and \code{incidence_input},
#' use the same time step as the zero-th time step.
#' @param interval_length Use with \code{estimation_method = "EpiEstim piecewise constant"}
#' Positive integer value.
#' Re is assumed constant over steps of size \code{interval_length}.
#'
#' @inheritParams inner_module
#' @inherit EpiEstim_wrapper
#'
.estimate_Re_EpiEstim_piecewise_constant <- function(incidence_input,
                                                     import_incidence_input = NULL,
                                                     minimum_cumul_incidence = 12,
                                                     interval_ends = NULL,
                                                     interval_length = 7,
                                                     mean_serial_interval = 4.8,
                                                     std_serial_interval = 2.3,
                                                     mean_Re_prior = 1,
                                                     output_HPD = FALSE) {

  .are_valid_argument_values(list(
    list(incidence_input, "module_input"),
    list(minimum_cumul_incidence, "non_negative_number"),
    list(interval_length, "positive_integer"),
    list(mean_serial_interval, "number"),
    list(std_serial_interval, "non_negative_number"),
    list(mean_Re_prior, "number")
  ))

  incidence_vector <- .get_values(incidence_input)

  if (sum(incidence_vector) < minimum_cumul_incidence) {
    stop("minimum_cumul_incidence parameter is set higher than total cumulative incidence.")
  }

  if(!is.null(import_incidence_input)) {
    .are_valid_argument_values(list(
      list(import_incidence_input, "module_input")
    ))

    incidence <- merge_outputs(output_list = list(local = incidence_input,
                                                  imported = import_incidence_input),
                               include_index = FALSE) %>%
      tidyr::replace_na(list(local = 0, imported = 0))

    incidence_length <- nrow(incidence)
    input_offset <- min(.get_offset(incidence_input), .get_offset(import_incidence_input))

  } else {
    incidence <- incidence_vector
    incidence_length <- length(incidence)
    input_offset <- .get_offset(incidence_input)
  }

  offset <- which(cumsum(incidence_vector) >= minimum_cumul_incidence)[1]
  # offset needs to be at least two for EpiEstim
  offset <- max(2, offset)
  right_bound <- incidence_length

  if (!is.null(interval_ends)) {
    .are_valid_argument_values(list(
      list(interval_ends, "integer_vector")
      ))

    # we make these be relative to input_offset
    interval_ends <- sort(interval_ends - input_offset)

    interval_ends <- interval_ends[interval_ends > offset & interval_ends <= right_bound]
  } else {
    interval_ends <- seq(from = offset + interval_length - 1, to = right_bound, by = interval_length)
    if (max(interval_ends) < right_bound) {
      interval_ends <- c(interval_ends, right_bound)
    }
  }

  if (length(interval_ends) < 1) {
    stop("No valid interval to estimate Re on.
         Check 'minimum_cumul_incidence', 'interval_ends' or 'interval_length' parameters.")
  }

  interval_starts <- c(offset, interval_ends[-length(interval_ends)] + 1)

  R_instantaneous <- suppressWarnings(EpiEstim::estimate_R(
    incid = incidence,
    method = "parametric_si",
    config = EpiEstim::make_config(
      list(
        mean_si = mean_serial_interval,
        std_si = std_serial_interval,
        t_start = interval_starts,
        t_end = interval_ends,
        mean_prior = mean_Re_prior
      )
    )
  )
  )

  additional_offset <- interval_starts[1] - 1

  replicate_estimates_on_interval <- function(estimates, interval_starts, interval_ends) {
    replicated_estimates <- unlist(lapply(
      seq_along(interval_starts),
      function(x) {
        rep(estimates[x], interval_ends[x] - interval_starts[x] + 1)
      }
    ))
    return(replicated_estimates)
  }

  Re_estimate <- replicate_estimates_on_interval(R_instantaneous$R$`Mean(R)`, interval_starts, interval_ends)
  Re_estimate <- .get_module_output(Re_estimate, input_offset, additional_offset)

  if (output_HPD) {
    Re_highHPD <- replicate_estimates_on_interval(R_instantaneous$R$`Quantile.0.975(R)`, interval_starts, interval_ends)
    Re_highHPD <- .get_module_output(Re_highHPD, input_offset, additional_offset)

    Re_lowHPD <- replicate_estimates_on_interval(R_instantaneous$R$`Quantile.0.025(R)`, interval_starts, interval_ends)
    Re_lowHPD <- .get_module_output(Re_lowHPD, input_offset, additional_offset)

    return(list(
      Re_estimate = Re_estimate,
      Re_highHPD = Re_highHPD,
      Re_lowHPD = Re_lowHPD
    ))
  } else {
    return(Re_estimate)
  }
}
