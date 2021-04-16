#TODO make class for input objects (list with values and index_offset fields)

#TODO make a utils that does the input checking (for before transformation and for when it's supposed to have been done)

#TODO make a utils function that converts individual module object into a tibble with a date

#TODO clarify the language between module input object and module output object

#TODO reorganize this file (maybe rename)

# Useful operator
`%!in%` <- Negate(`%in%`)

#' Merge multiple module outputs into tibble
#'
#' Output tibble from list of unsynced module outputs, with an optional date reference.
#' The optional reference date is the starting date of an input with offset 0.
#' In general, this will be the date corresponding to the first entry in the original incidence data.
#' If a reference date is provided, a date column is appended to the tibble,
#' with sequential dates generated with the time step specified by the optional time_step parameter.
#'
#'
#' @param output_list named list of module outputs
#' @param ref_date Date. Reference date.
#' @param time_step string. "day", "2 days", "week", "year"... (see seq.Date help page for details)
#'
#' @return tibble
#' @export
#'
#' @examples
#' #TODO add examples
merge_outputs <- function(output_list, ref_date = NULL, time_step = "day"){

  tibble_list <- lapply(1:length(output_list),
                        function(i) {
                          .make_tibble_from_output(output = output_list[[i]],
                                                  output_name = names(output_list)[i])
                        })

  merged_outputs <- plyr::join_all(tibble_list, by='index', type='full') %>%
                      dplyr::arrange(.data$index)

  if( !is.null(ref_date) ) {
    dates <- seq.Date(from = ref_date + min(merged_outputs$index), by = time_step, along.with = merged_outputs$index)
    merged_outputs$date <- dates
    merged_outputs <- dplyr::select(merged_outputs, date, tidyselect::everything())
  }

  merged_outputs <- dplyr::select(merged_outputs, -.data$index)

  return(merged_outputs)
}



#' Transform input data into a module input object
#'
#' #TODO specify what format data can take
#' @param data
#'
#' @return module input object
.get_module_input <- function(data) {
  #TODO properly check input format
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
.make_empty_module_output <- function(){
  return(list("values" = NA_real_, "index_offset" = 0))
}


#' Get values from module object
#'
#' This function must be adapted if the module_input_object implementation changes.
#'
#' @param module_object
#'
#' @return vector containing values
.get_values <- function(module_object) {
  if(is.list(module_object)) {
  return(module_object$values)
  } else {
    return(module_object)
  }
}


#' Get offset from module object
#'
#' This function must be adapted if the module_input_object implementation changes.
#'
#' @param module_object
#'
#' @return numeric scalar. Offset
.get_offset <- function(module_object) {
  if(is.list(module_object)) {
    return(module_object$index_offset)
  } else {
    return(0)
  }
}


#' Transform a result output of a module into a module output
#'
#' Also takes the module input object to calculate the offset of the output object.
#' The new offset is simply (module input offset) + (shift applied during module operations)
#' The shift applied during module operations is passed as "offset" parameter.
#'
#' @param results numeric vector
#' @param input module input object
#' @param offset integer
#'
#' @return module output object
.get_module_output <- function(results, input, offset = 0) {

  if(length(results) == 0) {
    return(.make_empty_module_output())
  }

  new_offset <- input$index_offset + offset

  while(is.na(results[1])) {
    results <- results[-1]
    new_offset <- new_offset + 1

    if(length(results) == 0) {
      return(.make_empty_module_output())
    }
  }

  return(list("values" = results, "index_offset" = new_offset))
}

#' Simplify output object if possible
#'
#' If offset is 0, return only vector containing values.
#'
#' @param output Module output object
#'
#' @return numeric vector or module output object. simplified output
.simplify_output <- function(output){
  if(.get_offset(output) == 0) {
    return(.get_values(output))
  } else {
    return(output)
  }
}


#' Get length of values vector in a module input object.
#'
#' @param input module input object.
#'
#' @return integer. length of values vector.
.get_input_length <- function(input) {
  return(length(.get_values(input)))
}


#' Convert module output object into tibble
#'
#' @param output module output object
#' @param output_name string. Name to be given to the values column
#'
#' @return tibble
.make_tibble_from_output <- function(output, output_name){

  tmp_output <- .get_module_input(output)
  indices <- seq(from = .get_offset(tmp_output), by = 1, length.out = length(.get_values(tmp_output)))

  return(dplyr::tibble(index = indices, !!output_name := .get_values(tmp_output)))
}

#TODO fill in documentation
#TODO allow for other types of time steps than days
#TODO stop exporting after removed from vignette
#' Generate delay data.
#'
#' This utility can be used to build toy examples to test functions dealing with empirical delay data.
#'
#' @param origin_date
#' @param n_time_steps
#' @param delay_ratio_start_to_end
#' @param time_step
#' @param distribution_initial_delay
#' @param seed
#'
#' @return data.frame. Simulated delay data.
#' @export
#'
#' @examples
#' #TODO add examples
generate_delay_data <- function(origin_date = as.Date("2020-02-01"),
                                 n_time_steps = 100,
                                 time_step = "day",
                                 delay_ratio_start_to_end = 2,
                                 distribution_initial_delay = list(name = "gamma", shape = 6, scale = 5),
                                 seed = NULL){

  if(!is.null(seed)) {
    set.seed(seed)
  }

  random_pop_size_walk <- cumsum(sample(c(-1, 0, 1), n_time_steps, TRUE))
  pop_size <- random_pop_size_walk - min(0, min(random_pop_size_walk)) + 1

  event_dates <- seq.Date(from = origin_date, length.out = n_time_steps, by = time_step)

  delays <- lapply(1:n_time_steps, function(i) {
    # Sample a number of draws from the gamma delay distribution, based on the pop size at time step i
    raw_sampled_delays <- .sample_from_distribution(distribution = distribution_initial_delay,
                                                    n = pop_size[i])
    # Multiply these samples by a factor that accounts for the linear inflation or deflation of delays.
    sampled_delays <- round(raw_sampled_delays * (1 + i/n_time_steps * (delay_ratio_start_to_end - 1)))

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
.print_vector <- function(a, digits = 3) {
  cat(paste0("c(",paste(round(a, digits = digits), collapse=","), ")"))
}



