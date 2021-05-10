#TODO make class for input objects (list with values and index_offset fields)
#TODO make a utils that does the input checking (for before transformation and for when it's supposed to have been done)
#TODO make a utils function that converts individual module object into a tibble with a date
#TODO clarify the language between module input object and module output object
#TODO reorganize this file (maybe rename)
#TODO write function that does inner addition of two incidence module inputs.

# Useful operator
`%!in%` <- Negate(`%in%`)

#List containing predefined accepted string inputs for exported functions, for parameters for which validity is tested using the.is_value_in_accepted_values_vector() function
accepted_parameter_value <- list(smoothing_method = c("LOESS"),
                                 deconvolution_method = c("Richardson-Lucy delay distribution"),
                                 estimation_method = c("EpiEstim sliding window"),
                                 bootstrapping_method = c("non-parametric block boostrap"),
                                 uncertainty_summary_method = c("original estimate - CI from bootstrap estimates", "bagged mean - CI from bootstrap estimates"))


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
                          index_col = "idx"){

  tibble_list <- lapply(1:length(output_list),
                        function(i) {
                          .make_tibble_from_output(output = output_list[[i]],
                                                   output_name = names(output_list)[i],
                                                   index_col = index_col)
                        })

  merged_outputs <- plyr::join_all(tibble_list, by= index_col, type='full') %>%
    dplyr::arrange(.data[[index_col]])

  if( !is.null(ref_date) ) {
    dates <- seq.Date(from = ref_date + min(merged_outputs[[index_col]]),
                      by = time_step,
                      along.with = merged_outputs[[index_col]])
    merged_outputs$date <- dates
    merged_outputs <- dplyr::select(merged_outputs, date, tidyselect::everything())
  }

  if(!include_index) {
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
                                     index_col = "idx"){

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

  dates <- seq.Date(from = ref_date + min(estimates[[index_col]]),
                    by = time_step,
                    along.with = estimates[[index_col]])
  estimates$date <- dates
  estimates <- dplyr::select(estimates, date, tidyselect::everything())

  if(!keep_index_col) {
    estimates <-  estimates %>%
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
#' @param module_object module object.
#'
#' @return vector containing \code{values} only
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
#' @param module_object module object.
#'
#' @return numeric scalar. \code{index_offset} of the \code{module_object}
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
#' The shift applied during module operations is the \code{offset} argument.
#'
#' @param results numeric vector containing output of module operations.
#' @param input module input object. Input originally given to module.
#' @param offset integer. Shift resulting from operations performed in module.
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

#TODO polish doc
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
correct_for_partially_observed_data <- function( incidence_data,
                                                 delay_distribution_final_report,
                                                 cutoff_observation_probability = 0.1 ) {

  #TODO validate cutoff_observation_probability argument
  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(delay_distribution_final_report, "delay_object", .get_input_length(incidence_data))))

  input <- .get_module_input(incidence_data)
  incidence_vector <- .get_values(input)

  #TODO clear up the mess: decide what kind of input delay_distribution_final_report can be and handle every case allowed
  delay_distribution_final_report_vector <- .get_delay_distribution(delay_distribution_final_report)

  #TODO build matrix beforehand (is not necessarily vector or matrix)
  if(NCOL(delay_distribution_final_report) == 1) {
    delay_distribution_matrix_final_report <- .get_matrix_from_single_delay_distr(delay_distribution_final_report_vector,
                                                                                  N=length(incidence_vector))
  } else {
    delay_distribution_matrix_final_report <- delay_distribution_final_report
  }

  Q_vector_observation_to_final_report <- apply(delay_distribution_matrix_final_report, MARGIN = 2, sum)

  #TODO improve this error
  if(any(is.na(Q_vector_observation_to_final_report)) || isTRUE(any(Q_vector_observation_to_final_report == 0, na.rm = FALSE))) {
    warning("Invalid delay_distribution_final_report argument.")
  }
  #TODO need to make sure that the matrix is the same size (as opposed to having extra columns leading)
  incidence_vector <- incidence_vector / Q_vector_observation_to_final_report

  # Now we cut off values at the end of the time series,
  # those dates for which the probability of having observed an event that happened on that date is too low
  # We define 'too low' as being below a 'cutoff_observation_probability'
  tail_values_below_cutoff <- which(rev(Q_vector_observation_to_final_report) < cutoff_observation_probability )

  if(length(tail_values_below_cutoff) == 0) {
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
    sampled_delays <- round(raw_sampled_delays * (1 + i/n_time_steps * (ratio_delay_end_to_start - 1)))

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

#' Return dots arguments as list.
#'
#' If there are no dots arguments, then return an empty list.
#'
#' @param ... dots arguments.
#'
#' @return list containing dots arguments or empty list
.get_dots_as_list <- function(...){
  if(...length() > 0) {
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
.get_arg_names <- function(func){
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
.get_shared_args <- function(func_list, dots_args){
  if(is.function(func_list)) {
    func_arg_names <- .get_arg_names(func_list)
  } else if(is.list(func_list)) {
    func_arg_names <- unlist(lapply(func_list, .get_arg_names))
  }
  return(dots_args[names(dots_args) %in% func_arg_names])
}


#' @description Utility function that checks if a specific user given parameter value is among the accepted ones, in which case it returns TRUE
#' Throws an error otherwise.
#'
#' @inherit validation_utility_params
.is_value_in_accepted_values_vector <- function(string_user_input, parameter_name){
  if(!is.character(string_user_input)){
    stop(paste("Expected parameter", parameter_name, "to be a string."))
  }
  if(!(string_user_input %in% accepted_parameter_value[[parameter_name]])){
    stop(paste("Expected parameter", parameter_name, "to have one of the following values:", toString(accepted_parameter_value[[parameter_name]]),". Given input was:", string_user_input))
  }
  return(TRUE)
}


#' @description Utility function that checks if a specific user given parameter value an accepted time_step, in which case it returns TRUE
#' An accepted time_step is considered to be:
#' <<A character string, containing one of "day", "week", "month", "quarter" or "year".
#' This can optionally be preceded by a (positive or negative) integer and a space, or followed by "s".>>
#' (from \url{https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/seq.Date})
#' @inherit validation_utility_params
#'
.is_value_valid_time_step <- function(string_user_input, parameter_name){
  if(!is.character(string_user_input)){
    stop(paste("Expected parameter", parameter_name, "to be a string."))
  }
  is_valid_time_step <- grepl("^([-+]?\\d+ )?(day|week|month|quarter|year)s?$", string_user_input)
  if(!is_valid_time_step){
    stop(paste("Expected parameter", parameter_name, "to be a character string, containing one of \"day\", \"week\", \"month\", \"quarter\" or \"year\". This can optionally be preceded by a (positive or negative) integer and a space, or followed by \"s\"."))
  }
  return(TRUE)
}


#' @description Utility function to determine whether an object is a numeric vector with all positive (or zero) values.
#'
#' @inherit validation_utility_params
#' @param vector vector to be tested
#'
#' @return boolean. TRUE if vector is a positive numeric vector. FALSE otherwise
.is_positive_numeric_vector <- function(vector){
  if(!is.vector(vector, mode="numeric")){
    return(FALSE)
  }
  if(!all(vector >= 0)){
    return(FALSE)
  }
  return(TRUE)
}


#' @description Utility function that checks if a user input is one of:
#' \itemize{
#'     \item a numeric vector with values > 0
#'     \item a list with two elements: \code{values} (a numeric vector with values > 0) and \code{index_offset} (an integer)
#' }
#' @inherit validation_utility_params
#' @param module_input_object the vector/list the user passed as a parameter, to be tested
#'
.is_valid_module_input <- function(module_input_object, parameter_name){
  if(is.list(module_input_object)){
    if("values" %!in% names(module_input_object)){
      stop(paste("When passed as a list,", parameter_name, "has to contain a $values element."))
    }

    if("index_offset" %!in% names(module_input_object)){
      stop(paste("When passed as a list,", parameter_name, "has to contain a $index_offset element."))
    }

    if(!.is_positive_numeric_vector(module_input_object$values)){
      stop(paste("The $values element of", parameter_name, "has to be a numeric vector with values greater or equal to 0."))
    }

    if(module_input_object$index_offset != as.integer(module_input_object$index_offset)){ #if index_offset is not an integer
      stop(paste("The $index_offset element of", parameter_name, "has to be an integer."))
    }

  } else if(is.numeric(module_input_object)){
    if(!.is_positive_numeric_vector(module_input_object)){
      stop(paste(parameter_name, "has to be a numeric vector with values greater or equal to 0."))
    }

  } else {
    stop(paste(parameter_name, "has to be either a numeric vector or a list."))
  }
  return(TRUE)
}

#' @description Utility function that checks if a given matrix is a valid delay distribution matrix.
#' For this, the matrix needs to fulfill the following conditions:
#' \itemize{
#'     \item is a numeric matrix
#'     \item has no values < 0
#'     \item is a lower triangular matrix
#'     \item no column sums up to more than 1
#'     \item no NA values
#'     \item the size of the matrix is greater than the length of the incidence data
#' }
#'
#' @inherit validation_utility_params
#' @param delay_matrix A matrix to be tested
#'
.check_is_delay_distribution_matrix <- function(delay_matrix, incidence_data_length){
  if(!is.matrix(delay_matrix) || !is.numeric(delay_matrix)){
    stop("The delay distribution object needs to be a numeric matrix.")
  }

  if(any(is.na(delay_matrix))){
    stop("The delay distribution matrix cannot contain any NA values.")
  }

  if(!all(delay_matrix >= 0)){
    stop("The delay distribution matrix needs to contain non-negative values.")
  }

  if(ncol(delay_matrix) != nrow(delay_matrix)){
    stop("The delay distribution matrix needs to be square.")
  }

  if(!all(delay_matrix == delay_matrix*lower.tri(delay_matrix, diag = TRUE))){ #check if matrix is lower triangular
    stop("The delay distribution matrix needs to be lower triangular.")
  }

  if(!all(colSums(delay_matrix) < 1)){
    stop("The delay distribution matrix is not valid. At least one column sums up to a value greater than 1.")
  }

  if(ncol(delay_matrix) < incidence_data_length){
    stop("The delay distribution matrix needs to have a greater size than the length of the incidence data.")
  }

  return(TRUE)

}

#' @description Utility function that checks whether a user input is a valid delay object. This means it can be one of the following:
#'      \itemize{
#'         \item a probability distribution vector: a numeric vector with no \code{NA} or negative values, whose entries sum up to 1
#'         \item an empirical delay data: a data frame with two columns: \code{event_date} and \code{report_delay}. The columns cannot contain \code{NA} values. \code{report_delay} only contains non-negative values
#'         \item a delay distribution matrix (as described in \code{\link{.check_is_delay_distribution_matrix}})
#'         \item a distribution object (e.g. list(name = 'gamma', scale = X, shape = Y))
#'      }
#' @inherit validation_utility_params
#' @param delay_object user inputted object to be tested
#'
.is_valid_delay_object <- function(delay_object, parameter_name, incidence_data_length){

  if(.is_numeric_vector(delay_object)){

    .check_is_probability_distr_vector(delay_object)

  } else if(is.data.frame(delay_object)){

    .check_is_empirical_delay_data(delay_object)

  } else if(is.matrix(delay_object)){

    .check_is_delay_distribution_matrix(delay_object, incidence_data_length)

  } else if(is.list(delay_object)){

    .is_valid_distribution(delay_object)

  } else {
    stop(paste("Invalid", parameter_name, "input.", parameter_name, "must be either:
         a numeric vector representing a discretized probability distribution,
         or a matrix representing discretized probability distributions,
         or a distribution object (e.g. list(name = 'gamma', scale = X, shape = Y)),
         or empirical delay data."))
  }
  return(TRUE)
}

#' @description  Utility function to check whether an object belongs to a particular class.
#' Wrapper function over \code{\link{.check_class}} needed because, being called from \code{\link{.are_valid_argument_values}},
#' the parameter name will not be the same as the one from the original function.
#'
#' @inherit validation_utility_params
#' @inherit .check_class
#'
.check_class_parameter_name <- function(object, proper_class, parameter_name, mode = "any"){
  tryCatch(
    {
      if(is.na(object)){
        stop("Object was NA") # This error message is never shown. Overwritten below.
      }
      .check_class(object, proper_class, mode)
    },
    error=function(error) {
      stop(paste("Expected parameter", parameter_name, "to be of type", proper_class, "and not NA."))
    }
  )
  return(TRUE)
}

#' @description Utility function to check whether an object is null or belongs to a particular class.
#'
#' @inherit validation_utility_params
#' @inherit .check_class
#'
.check_if_null_or_belongs_to_class <- function(object, proper_class, parameter_name, mode="any"){
  if(!is.null(object)){
    .check_class_parameter_name(object, proper_class, parameter_name, mode)
  }
  return(TRUE)
}


#' @description Utility function to check whether an object is a number.
#'
#' @inherit validation_utility_params
#' @param number The value to be tested
#'
.check_if_number <- function(number, parameter_name){
  if(!is.numeric(number)){
    stop(paste(parameter_name, "is expected to be a number."))
  }
  if(length(number) > 1){
    stop(paste(parameter_name, "is expected to be a number."))
  }
  return(TRUE)
}


#' @description Utility function to check whether an object is a positive number or 0.
#'
#' @inherit validation_utility_params
#' @inherit  .check_if_number
#'
.check_if_non_negative_number <- function(number, parameter_name){
  .check_if_number(number, parameter_name)

  if(number < 0){
    stop(paste(parameter_name, "is expected to be positive."))
  }

  return(TRUE)
}

#' @description Utility function that checks that the values the user passed when calling a function are valid.
#'
#' @inherit validation_utility_params
#' @param user_inputs A list of lists with two elements: the first is the value of the parameter to be tested. The second is the expected type of that parameter.
#'
.are_valid_argument_values <- function(user_inputs){
  for(i in 1:length(user_inputs)){
    user_input <- user_inputs[[i]][[1]]
    input_type <- user_inputs[[i]][[2]]
    parameter_name <- deparse(substitute(user_inputs)[[i+1]][[2]])
    if(length(user_inputs[[i]]) > 2){
      additional_function_parameter <- user_inputs[[i]][[3]]
    }

    switch (input_type,
            "smoothing_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
            "deconvolution_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
            "estimation_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
            "uncertainty_summary_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
            "bootstrapping_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
            "time_step" = .is_value_valid_time_step(user_input, parameter_name),
            "module_input" = .is_valid_module_input(user_input, parameter_name),
            "boolean" = .check_class_parameter_name(user_input,"logical", parameter_name),
            "delay_object" = .is_valid_delay_object(user_input, parameter_name, additional_function_parameter),
            "number" = .check_if_number(user_input, parameter_name),
            "non_negative_number" = .check_if_non_negative_number(user_input, parameter_name),
            "null_or_date" = .check_if_null_or_belongs_to_class(user_input, "Date", parameter_name),
            stop(paste("Checking function for type", input_type, "not found."))
    )
  }
  return(TRUE)
}
