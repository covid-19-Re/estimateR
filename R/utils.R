# TODO make class for input objects (list with values and index_offset fields)
# TODO make a utils that does the input checking (for before transformation and for when it's supposed to have been done)
# TODO make a utils function that converts individual module object into a tibble with a date
# TODO clarify the language between module input object and module output object
# TODO reorganize this file (maybe rename)
# TODO write function that does inner addition of two incidence module inputs.

# Useful operator
`%!in%` <- Negate(`%in%`)



#' Merge multiple module outputs into tibble
#'
#' Output tibble from list of unsynced module outputs, with an optional date column.
#' The optional \code{ref_date} argument is the starting date of an input with offset 0.
#' In general, this will be the date corresponding to the first entry in the original incidence data.
#' If a reference date is provided with \code{ref_date}, a date column is appended to the tibble,
#' with sequential dates generated with the time step specified by the \code{time_step} parameter.
#'
#'
#' @param output_list named list of module output objects.
#' @param include_index boolean. Include an index column in output?
#' @param index_col string. If \code{include_index} is \code{TRUE},
#' an index column named \code{index_col} is added to the output.
#' @inherit dating
#'
#' @return tibble
#' @export
merge_outputs <- function(output_list,
                          ref_date = NULL,
                          time_step = "day",
                          include_index = is.null(ref_date),
                          index_col = "idx") {
  .are_valid_argument_values(list(
    list(ref_date, "null_or_date"),
    list(time_step, "time_step"),
    list(index_col, "string"),
    list(include_index, "boolean")
  ))
  for (i in 1:length(output_list)) {
    .are_valid_argument_values(list(list(output_list[[i]], "module_input")))
  }

  tibble_list <- lapply(
    1:length(output_list),
    function(i) {
      .make_tibble_from_output(
        output = output_list[[i]],
        output_name = names(output_list)[i],
        index_col = index_col
      )
    }
  )

  merged_outputs <- plyr::join_all(tibble_list, by = index_col, type = "full") %>%
    dplyr::arrange(.data[[index_col]])

  if (!is.null(ref_date)) {
    dates <- seq.Date(
      from = ref_date + min(merged_outputs[[index_col]]),
      by = time_step,
      along.with = merged_outputs[[index_col]]
    )
    merged_outputs$date <- dates
    merged_outputs <- dplyr::select(merged_outputs, date, tidyselect::everything())
  }

  if (!include_index) {
    merged_outputs <- dplyr::select(merged_outputs, -.data[[index_col]])
  }

  return(merged_outputs)
}

#' Convert module output object into tibble
#'
#' @param output module output object.
#' @param output_name string. Name to be given to the \code{values} column
#' @param index_col string. Name of the index column included in the output.
#'
#' @return tibble
.make_tibble_from_output <- function(output,
                                     output_name,
                                     index_col = "idx") {
  .are_valid_argument_values(list(
    list(output, "module_input"),
    list(output_name, "string"),
    list(index_col, "string")
  ))

  tmp_output <- .get_module_input(output)
  indices <- seq(from = .get_offset(tmp_output), by = 1, length.out = length(.get_values(tmp_output)))

  return(dplyr::tibble(!!index_col := indices, !!output_name := .get_values(tmp_output)))
}


#' Add dates column to dataframe.
#'
#' @param estimates dataframe. Estimates.
#' @param keep_index_col boolean. Keep index column in output?
#' @inherit dating
#' @inherit uncertainty
#'
#' @return estimates dataframe with dates column.
.add_date_column <- function(estimates,
                             ref_date,
                             time_step,
                             index_col = "idx",
                             keep_index_col = FALSE) {
  .are_valid_argument_values(list(
    list(estimates, "estimates", index_col),
    list(ref_date, "date"),
    list(time_step, "time_step"),
    list(index_col, "string"),
    list(keep_index_col, "boolean")
  ))

  dates <- seq.Date(
    from = ref_date + min(estimates[[index_col]]),
    by = time_step,
    along.with = estimates[[index_col]]
  )

  estimates <- estimates %>%
    dplyr::arrange(.data[[index_col]]) %>%
    dplyr::mutate(date = dates)

  estimates <- dplyr::select(estimates, date, tidyselect::everything())

  if (!keep_index_col) {
    estimates <- estimates %>%
      dplyr::select(-.data[[index_col]])
  }

  return(estimates)
}

#' Transform input data into a module input object
#'
#' The input can be a list containing a \code{values} element
#' and an \code{index_offset} element, potentially among others.
#' It can also be a vector containing numeric values.
#'
#' @param data TODO specify what format data can take
#'
#' @return module input object.
#' List with a \code{values} and an \code{index_offset} element.
#'
.get_module_input <- function(data) {
  .are_valid_argument_values(list(list(data, "module_input")))

  if (is.list(data)) {
    values <- as.double(data$values)
    index_offset <- data$index_offset
  } else {
    values <- as.double(data)
    index_offset <- 0
  }
  return(list("values" = values, "index_offset" = index_offset))
}

#' Make empty module output object
#'
#' @return empty module output object
.make_empty_module_output <- function() {
  return(list("values" = NA_real_, "index_offset" = 0))
}


#' Get values from module object
#'
#' This function must be adapted if the module_input_object implementation changes.
#'
#' @param module_object module object.
#'
#' @return vector containing \code{values} only
.get_values <- function(module_object) {
  .are_valid_argument_values(list(list(module_object, "module_input")))
  if (is.list(module_object)) {
    return(module_object$values)
  } else {
    return(module_object)
  }
}


#' Get offset from module object
#'
#' This function must be adapted if the module_input_object implementation changes.
#'
#' @param module_object module object.
#'
#' @return numeric scalar. \code{index_offset} of the \code{module_object}
.get_offset <- function(module_object) {
  .are_valid_argument_values(list(list(module_object, "module_input")))
  if (is.list(module_object)) {
    return(module_object$index_offset)
  } else {
    return(0)
  }
}


#' Transform a result output of a module into a module output
#'
#' Also takes the module input object to calculate the offset of the output object.
#' The new offset is simply (module input offset) + (shift applied during module operations)
#' The shift applied during module operations is the \code{offset} argument.
#'
#' @param results numeric vector containing output of module operations.
#' @param input module input object. Input originally given to module.
#' @param offset integer. Shift resulting from operations performed in module.
#'
#' @return module output object
.get_module_output <- function(results, input, offset = 0) {
  .are_valid_argument_values(list(
    list(results, "numeric_vector"),
    list(input, "module_input"),
    list(offset, "integer")
  ))
  if (length(results) == 0) {
    return(.make_empty_module_output())
  }

  new_offset <- input$index_offset + offset

  while (is.na(results[1])) {
    results <- results[-1]
    new_offset <- new_offset + 1

    if (length(results) == 0) {
      return(.make_empty_module_output())
    }
  }

  return(list("values" = results, "index_offset" = new_offset))
}

# TODO doc
#' Title
#'
#' @param input_a
#' @param input_b
#'
#' @return
#' @export
inner_addition <- function(input_a, input_b) {
  .are_valid_argument_values(list(
    list(input_a, "module_input"),
    list(input_b, "module_input")
  ))

  length_a <- .get_input_length(input_a)
  length_b <- .get_input_length(input_b)

  offset_a <- .get_offset(input_a)
  offset_b <- .get_offset(input_b)

  inner_offset <- max(offset_a, offset_b)
  length_addition <- min(length_a - (inner_offset - offset_a), length_b - (inner_offset - offset_b))

  inner_a <- .get_values(input_a)[seq(from = inner_offset - offset_a + 1, by = 1, length.out = length_addition)]
  inner_b <- .get_values(input_b)[seq(from = inner_offset - offset_b + 1, by = 1, length.out = length_addition)]

  return(.get_module_input(list(values = inner_a + inner_b, index_offset = inner_offset)))
}

# TODO doc
#' Title
#'
#' @param input_a
#' @param input_b
#'
#' @return
#' @export
left_addition <- function(input_a, input_b) {
  .are_valid_argument_values(list(
    list(input_a, "module_input"),
    list(input_b, "module_input")
  ))

  offset_a <- .get_offset(input_a)
  offset_b <- .get_offset(input_b)

  min_offset <- min(offset_a, offset_b)
  padded_input_a <- leftpad_input(input_a, min_offset, padding_value = 0)
  padded_input_b <- leftpad_input(input_b, min_offset, padding_value = 0)

  length_a <- .get_input_length(padded_input_a)
  length_b <- .get_input_length(padded_input_b)

  length_addition <- min(length_a, length_b)

  values_a <- .get_values(padded_input_a)[1:length_addition]
  values_b <- .get_values(padded_input_b)[1:length_addition]

  return(.get_module_input(list(values = values_a + values_b, index_offset = min_offset)))
}

# TODO doc
leftpad_input <- function(input, new_offset, padding_value = 0) {
  .are_valid_argument_values(list(
    list(input, "module_input"),
    list(new_offset, "number"),
    list(padding_value, "number")
  ))

  if (new_offset >= .get_offset(input)) {
    return(input)
  } else {
    padded_values <- c(rep(padding_value, length.out = .get_offset(input) - new_offset), .get_values(input))
    return(list(values = padded_values, index_offset = new_offset))
  }
}

# TODO test with matrix delay
# TODO polish doc
# TODO redoc
#' Correct incidence data for yet-to-be-observed fraction of events
#'
#' Use this function to correct the tail of an incidence timeseries
#' if incidence was collected following a subsequent observation event.
#' For instance, if the incidence represents people starting to show symptoms of a disease
#' (dates of onset of symptoms), the data would typically have been collected among
#' individuals whose case was confirmed via a test.
#' If so, among all events of onset of symptoms, only those who had time to be
#' confirmed by a test were reported.
#' Thus, close to the present, there is an underreporting of onset of symptoms events.
#' In order to account for this effect, this function divides each incidence value
#' by the probability of an event happening at a particular timestep to have been observed.
#' Typically, this correction only affects the few most recent datapoints.
#'
#' @param delay_distribution_final_report TODO refactor to estimateR
#' Distribution of the delay between the events collected in the incidence data
#' and the a posteriori observations of these events.
#' @param cutoff_observation_probability value between 0 and 1.
#' Only datapoints for timesteps that have a probability to be observed higher
#' than \code{cutoff_observation_probability} are kept.
#' The few datapoints with a lower probability to be observed are trimmed off
#' the tail of the timeseries.
#' @inheritParams module_structure
#'
#' @return module output object
#' @export
correct_for_partially_observed_data <- function(incidence_data,
                                                delay_distribution_final_report,
                                                cutoff_observation_probability = 0.1, # TODO set to a higher value i.e. 0.25 by default
                                                ref_date = NULL,
                                                time_step = "day",
                                                ...) {
  .are_valid_argument_values(list(
    list(incidence_data, "module_input"),
    list(delay_distribution_final_report, "delay_object", .get_input_length(incidence_data)),
    list(cutoff_observation_probability, "numeric_between_zero_one"),
    list(ref_date, "null_or_date"),
    list(time_step, "time_step")
  ))

  input <- .get_module_input(incidence_data)
  incidence_vector <- .get_values(input)

  dots_args <- .get_dots_as_list(...)

  delay_distribution_final_report <- do.call(
    ".get_delay_distribution",
    c(
      list(
        delay = delay_distribution_final_report,
        n_report_time_steps = length(incidence_vector),
        ref_date = ref_date,
        time_step = time_step
      ),
      .get_shared_args(
        list(
          .get_delay_distribution,
          get_matrix_from_empirical_delay_distr,
          build_delay_distribution
        ),
        dots_args
      )
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
    warning("Invalid delay_distribution_final_report argument.")
  }
  # TODO need to make sure that the matrix is the same size (as opposed to having extra columns leading)
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

#' Simplify output object if possible
#'
#' If offset is 0, return only vector containing values.
#' If offset is not zero then \code{output} is returned as is.
#'
#' @param output Module output object
#'
#' @return numeric vector or module output object.
.simplify_output <- function(output) {
  .are_valid_argument_values(list(list(output, "module_input")))


  if (.get_offset(output) == 0) {
    return(.get_values(output))
  } else {
    return(output)
  }
}

.simplify <- function(output_list,
                      ref_date,
                      time_step = "day") {
  merged_outputs <- merge_outputs(
    output_list = output_list,
    ref_date = ref_date,
    time_step = time_step
  )
}

#' Get length of values vector in a module input object.
#'
#' @param input module input object.
#'
#' @return integer. length of values vector.
.get_input_length <- function(input) {
  .are_valid_argument_values(list(list(input, "module_input")))
  return(length(.get_values(input)))
}

#' Generate delay data.
#'
#' This utility can be used to build toy examples to test functions dealing with empirical delay data.
#' It is very basic in what it simulates.
#' A random walk is simulated over \code{n_time_steps}, representing the incidence through time.
#' The result of this simulation is offset so that all values are positive.
#' Then, for each time step, \code{n} samples from a delay distribution are taken,
#' with \code{n} being the incidence value at this time step.
#' The random draws are then multiplied by a factor (>1 or <1) to simulate
#' a gradual shift in the delay distribution through time.
#' This multiplication factor is calculated
#' by linearly interpolating between 1 (at the first time step),
#' and \code{delay_ratio_start_to_end} linearly,
#' from 1 at the first time step to \code{ratio_delay_end_to_start}
#' at the last time step.
#'
#' @param origin_date Date of first infection.
#' @param n_time_steps interger. Number of time steps to generate delays over
#' @param ratio_delay_end_to_start numeric value.
#' Shift in delay distribution from start to end.
#' @param distribution_initial_delay Distribution in list format.
#' @param seed integer. Optional RNG seed.
#' @inherit dating
#'
#' @return dataframe. Simulated delay data.
.generate_delay_data <- function(origin_date = as.Date("2020-02-01"),
                                 n_time_steps = 100,
                                 time_step = "day",
                                 ratio_delay_end_to_start = 2,
                                 distribution_initial_delay = list(name = "gamma", shape = 6, scale = 5),
                                 seed = NULL) {
  .are_valid_argument_values(list(
    list(origin_date, "date"),
    list(n_time_steps, "positive_integer"),
    list(time_step, "time_step"),
    list(ratio_delay_end_to_start, "number"),
    list(distribution_initial_delay, "distribution"),
    list(seed, "null_or_int")
  ))

  if (!is.null(seed)) {
    set.seed(seed)
  }

  random_pop_size_walk <- cumsum(sample(c(-1, 0, 1), n_time_steps, TRUE))
  pop_size <- random_pop_size_walk - min(0, min(random_pop_size_walk)) + 1

  event_dates <- seq.Date(from = origin_date, length.out = n_time_steps, by = time_step)

  delays <- lapply(1:n_time_steps, function(i) {
    # Sample a number of draws from the gamma delay distribution, based on the pop size at time step i
    raw_sampled_delays <- .sample_from_distribution(
      distribution = distribution_initial_delay,
      n = pop_size[i]
    )
    # Multiply these samples by a factor that accounts for the linear inflation or deflation of delays.
    sampled_delays <- round(raw_sampled_delays * (1 + i / n_time_steps * (ratio_delay_end_to_start - 1)))

    return(tibble::tibble(event_date = event_dates[i], report_delay = sampled_delays))
  })

  return(dplyr::bind_rows(delays))
}

#' Utility function to print vectors in a copy-pastable format
#'
#' @param a numeric vector
#' @param digits integer. Number of digits to round.
#'
#' @return Print vector as string
.print_vector <- function(a, digits = 2) {
  cat(paste0("c(", paste(round(a, digits = digits), collapse = ","), ")"))
}

#' Return dots arguments as list.
#'
#' If there are no dots arguments, then return an empty list.
#'
#' @param ... dots arguments.
#'
#' @return list containing dots arguments or empty list
.get_dots_as_list <- function(...) {
  if (...length() > 0) {
    dots_args <- list(...)
  } else {
    dots_args <- list()
  }
  return(dots_args)
}

#' Get the names of the arguments of a function
#'
#' @param func function
#'
#' @return names of the arguments of \code{func}
.get_arg_names <- function(func) {
  return(names(formals(func)))
}

#' Get the arguments which apply to a function among a given list of arguments.
#'
#' This function is used to find which arguments should be passed
#' to a list of functions among the dot arguments passed to a higher-level function.
#'
#' @param func_list list of functions
#' @param dots_args list of arguments
#'
#' @return elements of \code{dots_args} which can be passed to \code{func_list}
.get_shared_args <- function(func_list, dots_args) {
  if (is.function(func_list)) {
    func_arg_names <- .get_arg_names(func_list)
  } else if (is.list(func_list)) {
    func_arg_names <- unlist(lapply(func_list, .get_arg_names))
  }
  return(dots_args[names(dots_args) %in% func_arg_names])
}
