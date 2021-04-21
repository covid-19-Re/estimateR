#TODO make class for input objects (list with values and index_offset fields)

#TODO make a utils that does the input checking (for before transformation and for when it's supposed to have been done)

#TODO make a utils function that converts individual module object into a tibble with a date

#TODO clarify the language between module input object and module output object

#TODO reorganize this file (maybe rename)

# Useful operator
`%!in%` <- Negate(`%in%`)

#List containing predefined accepted string inputs for exported functions, for parameters for which validity is tested using the.is_value_in_accepted_values_vector() function
accepted_parameter_value <- list(smoothing_method = c("LOESS"),
                                 deconvolution_method = c("Richardson-Lucy delay distribution"),
                                 estimation_method = c("EpiEstim sliding window"),
                                 bootstrapping_method = c("non-parametric block boostrap"))


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


#' Utility function that checks if a specific user given parameter value is among the accepted ones, in which case it returns TRUE
#' Throws an error otherwise.
#' @param string_user_input string containing the value that the user passed for the tested parameter
#' @param parameter_name string containing the name of the parameter to be tested
#'
#' @return a boolean value. (TRUE if string_user_input is an accepted value. Throws an error otherwise)
# TODO maybe add "" or , in between accepted parameter values in error message
.is_value_in_accepted_values_vector <- function(string_user_input, parameter_name){
  if(!is.character(string_user_input)){
    stop(paste("Expected parameter", parameter_name, "to be a string."))
  }
  if(!(string_user_input %in% accepted_parameter_value[[parameter_name]])){
    stop(paste("Expected parameter", parameter_name, "to have one of the following values:", toString(accepted_parameter_value[[parameter_name]]),"."))
  }
  return(TRUE)
}


#' Utility function that checks if a specific user given parameter value an accepted time_step, in which case it returns TRUE
#' Throws an error otherwise.
#' @param string_user_input string containing the value that the user passed for the tested parameter
#' @param parameter_name string containing the name of the parameter to be tested
#' 
#' @return a boolean value. (TRUE if string_user_input is an accepted time_step. Throws an error otherwise)
#' An accepted time_step is considered to be: <<A character string, containing one of "day", "week", "month", "quarter" or "year". This can optionally be preceded by a (positive or negative) integer and a space, or followed by "s".>> (from https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/seq.Date)
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


#' Utility function to determine whether an object is a numeric vector with all positive (or zero) values.
#'
#' @param vector vector to be tested
#'
#' @return boolean. TRUE if vector is a positive numeric vector. FALSE otherwise

.is_positive_numeric_vector <- function(vector){
  if(!is.numeric(vector)){
    return(FALSE)
  }
  if(!all(vector >= 0)){
    return(FALSE)
  }
  return(TRUE)
}


#' Utility function that checks if a user input respects one the following criteria:
#' - is either a numeric vector with values > 0
#' - or is a list with two elements: $values: numeric vector with values > 0
#'                                   $index_offset: integer 
#'
#' @param user_input the vector/list the user passed as a parameter, to be tested
#'
#' @return TRUE if user_input is a valid module input. Throws an informative error otherwise.
.is_valid_module_input <- function(module_object, parameter_name){
  if(is.list(module_object)){
    if("values" %!in% names(module_object)){
      stop(paste("When passed as a list,", parameter_name, "has to contain a $values element."))
    }
    
    if("index_offset" %!in% names(module_object)){
      stop(paste("When passed as a list,", parameter_name, "has to contain a $index_offset element."))
    } 
    
    if(!.is_positive_numeric_vector(module_object$values)){
      stop(paste("The $values element of", parameter_name, "has to be a numeric vector with values greater or equal to 0."))
    }
    
    if(module_object$index_offset != as.integer(module_object$index_offset)){ #if index_offset is not an integer
      stop(paste("The $index_offset element of", parameter_name, "has to be an integer."))
    } 
    
  } else if(is.numeric(module_object)){
    if(!.is_positive_numeric_vector(module_object)){
      stop(paste(parameter_name, "has to be a numeric vector with values greater or equal to 0."))
    }
    
  } else {
    stop(paste(parameter_name, "has to be either a numeric vector or a list."))
  }
  return(TRUE)
}

#' TODO: fill in 
#' TODO: can it have values on diagonal?
#' TODO: test
#' question: ask whether to keep throw error behavior or add a throw=True argument to the function?
#' @param delay_matrix 
#'
#' @return
.is_delay_distribution_matrix <- function(delay_matrix, incidence_data_length){
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
    
  if(!all(delay_matrix == delay_matrix*lower.tri(delay_matrix))){ #check if matrix is lower triangular
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

#' TODO fill in
#'
#' @param object 
#' @param proper_class 
#' @param parameter_name 
#'
#' @return
#' @export
#'
#' @examples
.check_class_parameter_name <- function(object, proper_class, parameter_name){
  tryCatch(
    {
      if(is.na(object)){
        stop("Object was NA") # This error message is never shown. Overwritten below. 
      }
      .check_class(object, proper_class)
    },
    error=function(error) {
      stop(paste("Expected parameter", parameter_name, "to be of type", proper_class, "and not NA."))
    }
  )    
  return(TRUE)
}

#' Utility function that checks that the values the user passed when calling a function are valid
#' Returns TRUE if all checks were passed.
#' @param user_inputs list of all arguments with which the tested function was called (can be obtain via "as.list(environment()")
#' 
#' @return TRUE if all checks were passed. Throws an error otherwise
.are_valid_argument_values <- function(user_inputs){
  for(i in 1:length(user_inputs)){
    user_input <- user_inputs[[i]]$user_input 
    input_type <- user_inputs[[i]]$input_type
    parameter_name <- user_inputs[[i]]$parameter_name
    switch (input_type,
        "smoothing_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
        "deconvolution_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
        "estimation_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
        "bootstrapping_method" = .is_value_in_accepted_values_vector(user_input, parameter_name),
        "time_step" = .is_value_valid_time_step(user_input, parameter_name),
        "module_input" = .is_valid_module_input(user_input, parameter_name),
        "boolean" = .check_class_parameter_name(user_input,"logical", parameter_name),
        "empirical_delay_data" = .check_is_empirical_delay_data(user_input)
    )
  }
  return(TRUE)
}

