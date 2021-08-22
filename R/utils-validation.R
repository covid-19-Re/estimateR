#List containing predefined accepted string inputs for exported functions, for parameters for which validity is tested using the.is_value_in_accepted_values_vector() function
accepted_parameter_value <- list(smoothing_method = c("LOESS"),
                                 deconvolution_method = c("Richardson-Lucy delay distribution"),
                                 estimation_method = c("EpiEstim sliding window", "EpiEstim piecewise constant"),
                                 bootstrapping_method = c("non-parametric block boostrap"),
                                 function_prefix = c("d", "q", "p", "r"),
                                 uncertainty_summary_method = c("original estimate - CI from bootstrap estimates", "bagged mean - CI from bootstrap estimates"),
                                 fit = c("none", "gamma"))

#' Check that an object represents a probability distribution.
#'
#' To pass the check:
#' 1) the object must be a numeric vector
#' 2) its elements must sum to 1
#' 3) it must not contain any strictly-negative value.
#' 4) (optionally) it must not contain NAs.
#'
#' @param distribution Input for which we need to check that it is a proper probability distribution.
#' @param tolerate_NAs Can the distribution contain NA values?
#' @param tolerance_on_sum Numeric tolerance in checking that vector elements sum to 1.
#'
#' @inherit validation_utility_params
.check_is_probability_distr_vector <- function(distribution, tolerate_NAs = FALSE, tolerance_on_sum = 1E-2, parameter_name = deparse(substitute(distribution))) {

  .check_class_parameter_name(distribution, "vector", parameter_name, mode = "numeric")
  if( !tolerate_NAs && any(is.na(distribution)) ) {
    stop("Not a proper delay distribution vector. Contains one or more NAs.")
  }

  if( !isTRUE(all.equal(1, sum(distribution, na.rm = TRUE ), tolerance = tolerance_on_sum))){
    stop("Not a proper delay distribution vector. Does not sum to 1.")
  }

  if( any(distribution  < 0, na.rm = TRUE)) {
    stop("Not a proper delay distribution vector. Contains negative values.")
  }

  return(TRUE)
}


#' Check whether the class of an object is as expected
#'
#' @param object An object whose class needs checking,
#' @param proper_class A string describing the desired class of \code{object}.
#' @param mode Optional. A string describing the desired mode of \code{object}.
#' Use only if \code{proper_class} is \code{vector}. Mode cannot be \code{Date}.
#' Use \code{proper_class = "Date"} for checking class of \code{Date vector}.
#'
#' @return TRUE if no error is thrown.
.check_class <- function(object, proper_class, mode = "any"){
  if ("character" %!in% class(proper_class) || length(proper_class) > 1 ) {
    stop("'proper_class' must be a single string.")
  }

  if ("character" %!in% class(mode) || length(mode) > 1 ) {
    stop("'mode' must be a single string.")
  }

  if(proper_class == "vector") {
    if(mode == "Date") {
      stop("Mode cannot be 'Date'.")
    }

    if(!is.vector(object, mode = mode)) {
      stop(paste0(deparse(substitute(object)), " must be a ", mode, " vector."))
    }

    return(TRUE)
  }

  # validation function
  is_proper_class <- get(paste0("is.", proper_class), envir = loadNamespace("lubridate")) # need lubridate in case proper_class is Date

  if (!is_proper_class(object)) {
    # deparse(substitute(...)) lets you do basically the reverse of get(..)
    stop(paste0(deparse(substitute(object)), " must be a ", proper_class, "."))
  }

  return(TRUE)
}

#TODO add checks for other distributions (lognormal, uniform, weibull, truncated_normal,...)
#TODO reconsider if can return FALSE
#' Check if valid distribution list
#'
#'
#'
#' @inheritParams distribution
#' @inheritParams validation_utility_params
#'
#'
#' @return boolean. Returns FALSE if parameter values return an improper distribution (if gamma distr). Throws an error if not a list, or not a list with the appropriate elements. Returns TRUE otherwise.
.is_valid_distribution <- function(distribution, parameter_name = deparse(substitute(distribution))){

  .check_class_parameter_name(distribution, "list", parameter_name)

  if (!"name" %in% names(distribution)) {
    stop("Missing distribution name. Include a 'name' element in distribution.")
  }

  distribution_name <- distribution[["name"]]

  density_function_name <- paste0("d", distribution_name)
  density_function <- try(get(density_function_name, envir = loadNamespace("stats")),
                          silent = TRUE)

  if (class(density_function) == "try-error") {
    stop(paste("The ", density_function_name, " function must be defined in the 'stats' package."))
  }

  distribution_parms <- .get_distribution_parms(distribution, density_function)

  if (length(distribution_parms) == 0) {
    stop("Missing distribution parameters.")
  }

  # Check if parameter values are pathological.
  if (distribution_name == "gamma") {

    if (distribution_parms[["shape"]] < 0) {
      return(FALSE)
    } else if ("scale" %in% names(distribution_parms) && distribution_parms[["scale"]] <= 0) {
      return(FALSE)
    } else {
      return(TRUE)
    }

  }

  return(TRUE)
}

#' Check if input is in the proper empirical delay data format
#'
#' If the \code{delay} input is not a dataframe, return \code{FALSE}.
#' Otherwise, an error is thrown if \code{delay} does not follow the expected format.
#'
#' @inherit empirical_delay_data_format details
#' @inherit validation_utility_params
#' @param delay object to be tested
#'
#' @return boolean. \code{TRUE} if the input is a dataframe in the proper format.
.check_is_empirical_delay_data <- function(delay, parameter_name = deparse(substitute(distribution))){
  if(is.data.frame(delay)) {

    if("event_date" %!in% colnames(delay)) {
      stop("Missing 'event_date' column in dataframe.")
    }
    .check_class_parameter_name(delay$event_date, "Date", parameter_name)

    if ("report_delay" %!in% colnames(delay)) {
      stop("Missing 'report_delay' column in dataframe.")
    }
    .check_class_parameter_name(delay$report_delay, "vector", parameter_name, mode = "numeric")

    if( any(is.na(delay$event_date)) || any(is.na(delay$report_delay)) ){
      stop("Empirical delay data contains NA values.")
    }

    if( any(delay$report_delay < 0) ){
      stop("'report_delay' column contains negative values.")
    }

    return(TRUE)

  } else {
  return(FALSE)
  }
}


#' Check if object is numeric vector.
#'
#' @param object Any object.
#'
#' @return TRUE if numeric vector, FALSE otherwise
.is_numeric_vector <- function(object){
  return(is.vector(object, mode = "numeric"))
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

#' Utility functions for input validity.
#' 
#' @description Utility function to determine whether an object is a numeric vector with all positive (or zero) values.
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

.is_list_of_outputs <- function(output_list){
  if(!is.list(output_list)) {
    return(FALSE)
  }

  check_if_simple_output <- try(.is_valid_module_input(output_list, deparse(substitute(output_list))),
                                silent = TRUE)

  # Return FALSE if input is an output object itself
  if(!("try-error" %in% class(check_if_simple_output))) {
    return(FALSE)
  }

  for(i in 1:length(output_list)) {
    test_output_i <- try(.is_valid_module_input(output_list[[i]], names(output_list)[i]),
                         silent = TRUE)
    if("try-error" %in% class(test_output_i)) {
      return(FALSE)
    }
  }

  return(TRUE)
}



#TODO reconsider whether we need the incidence_data_length here.
# Is it acceptable if dim(matrix) > incidence data length?
# And is it needed to check whether ncol(delay_matrix) < incidence_data_length
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
.check_is_delay_distribution_matrix <- function(delay_matrix, incidence_data_length, parameter_name){
  if(!is.matrix(delay_matrix) || !is.numeric(delay_matrix)){
    stop(paste(parameter_name, "needs to be a numeric matrix."))
  }

  if(any(is.na(delay_matrix))){
    stop(paste(parameter_name, "cannot contain any NA values."))
  }

  if(!all(delay_matrix >= 0)){
    stop(paste(parameter_name, "needs to contain non-negative values."))
  }

  if(ncol(delay_matrix) != nrow(delay_matrix)){
    stop(paste(parameter_name, "needs to be a square matrix."))
  }

  if(!all(delay_matrix == delay_matrix*lower.tri(delay_matrix, diag = TRUE))){ #check if matrix is lower triangular
    stop(paste(parameter_name, "needs to be a lower triangular matrix."))
  }

  if(!all(colSums(delay_matrix) <= 1)){
    stop(paste(parameter_name, "is not a valid delay distribution matrix. At least one column sums up to a value greater than 1."))
  }

  if(ncol(delay_matrix) < incidence_data_length){
    stop(paste(parameter_name,"needs to have a greater size than the length of the incidence data."))
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

    .check_is_probability_distr_vector(delay_object, parameter_name = parameter_name)

  } else if(is.data.frame(delay_object)){

    .check_is_empirical_delay_data(delay_object, parameter_name)

  } else if(is.matrix(delay_object)){

    .check_is_delay_distribution_matrix(delay_object, incidence_data_length, parameter_name)

  } else if(is.list(delay_object)){

    .is_valid_distribution(delay_object, parameter_name)

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
      if(length(object) == 1 && is.na(object)){
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

#' @description Utility function to check whether an object is a strictly positive number
#'
#' @inherit validation_utility_params
#' @inherit  .check_if_number
#'
.check_if_positive_number <- function(number, parameter_name){
  .check_if_number(number, parameter_name)

  if(number <= 0){
    stop(paste(parameter_name, "is expected to be strictly positive."))
  }

  return(TRUE)

}

#' @description Utility function to check whether an object is an integer
#'
#' @inherit validation_utility_params
#' @param number The value to be tested
#'
.check_if_integer <- function(number, parameter_name){
  if(as.integer(number) != number){ # did not use .check_class_parameter_name since is.integer(1) returns false
    stop(paste(parameter_name, "needs to be an integer."))
  }
  return(TRUE)
}

#' @description Utility function to check whether an object is an integer or null
#'
#' @inherit validation_utility_params
#' @inherit  .check_if_integer
#'
.check_if_null_or_integer <- function(number, parameter_name){
  if(!is.null(number)){
    .check_if_integer(number, parameter_name)
  }
  return(TRUE)
}


#' @description Utility function to check whether an object is a strictly positive integer
#'
#' @inherit validation_utility_params
#' @inherit  .check_if_integer
#'
.check_if_positive_integer <- function(number, parameter_name){
  .check_if_positive_number(number, parameter_name)
  .check_if_integer(number, parameter_name)
}

#' @description Utility function to check whether an object is a number that belongs to a given interval
#'
#' @inherit validation_utility_params
#' @inherit  .check_if_number
#'
.check_is_numeric_in_interval  <- function(user_input, parameter_name, interval_start, interval_end){
  .check_if_number(user_input, parameter_name)
  if(user_input < interval_start || user_input > interval_end){
    stop(paste0("Expected ", parameter_name, " to be in interval [", interval_start,", ", interval_end, "]."))
  }
  return(TRUE)
}


#' @description Utility function to check whether an object is a valid estimates object.
#' It must be a dataframe and have an index column named \code{index_col_name} that doesn't contain any \code{NA} values.
#' @inherit validation_utility_params
#' @param index_col_name string. Name of the index column in the \code{user_input} dataframe.
#'
.check_is_estimate  <- function(user_input, parameter_name, index_col_name){
  .check_class_parameter_name(user_input, "data.frame", parameter_name)
 
  if(index_col_name %!in% names(user_input)){
    stop(paste("Missing index column. No column named ", index_col_name, "in", parameter_name))
  }

  if(any(is.na(user_input[[index_col_name]]))) {
    stop(paste("NA value(s) in column", index_col_name, "in", parameter_name))
  }
  return(TRUE)
}


#' @description Utility function to check whether an object a valid bootstrap estimates object. 
#' It has to be a valid estimates object, and to have the columns specified by \code{col_names}.
#' 
#' @inherit validation_utility_params
#' @param col_names vector. Contains the column names of \code{index_col}, \code{bootstrap_id_col} and \code{Re_estimate_col}, as described by the \code{summarise_uncertainty} function.
#'
.check_is_bootstrap_estimate  <- function(user_input, parameter_name, col_names){
  Re_estimate_col = col_names[1]
  bootstrap_id_col = col_names[2]
  index_col = col_names[3]

  .check_is_estimate(user_input, parameter_name, index_col)

  for(i in 1:2){#the bootstrap_id column name and Re_estimate column name; index column is already checked by .check_is_estimate
    if(col_names[i] %!in% names(user_input)) {
      stop(paste0("Missing ", col_names[i], " column in 'bootstrapped estimates' argument,
                or '", col_names[i],"' was not set to the corresponding column name."))
    }
  }

  return(TRUE)
}

#' Utility functions for input validity.
#' 
#' @description Utility function that checks that the values the user passed when calling a function are valid.
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
            "null_or_int" = .check_if_null_or_integer(user_input, parameter_name),
            "positive_integer" = .check_if_positive_integer(user_input, parameter_name),
            "positive_number" = .check_if_positive_number(user_input, parameter_name),
            "string" = .check_if_null_or_belongs_to_class(user_input, "character", parameter_name),
            "date" = .check_class_parameter_name(user_input, "Date", parameter_name),
            "integer" = .check_if_integer(user_input, parameter_name),
            "distribution" = .is_valid_distribution(user_input, parameter_name),
            "numeric_between_zero_one" = .check_is_numeric_in_interval(user_input, parameter_name, 0, 1),
            "function_prefix" = .is_value_in_accepted_values_vector(user_input, parameter_name),
            "numeric_vector" = .check_class_parameter_name(user_input, "vector", parameter_name, "numeric"),
            "probability_distr_vector" = .check_is_probability_distr_vector(user_input, parameter_name = parameter_name),
            "probability_distr_vector_high_tolerance" = .check_is_probability_distr_vector(user_input, parameter_name = parameter_name, tolerance_on_sum = 1E-2),
            "probability_distr_matrix" = .check_is_delay_distribution_matrix(user_input, additional_function_parameter, parameter_name),
            "empirical_delay_data" = .check_is_empirical_delay_data(user_input, parameter_name),
            "estimates" = .check_is_estimate(user_input, parameter_name, additional_function_parameter),
            "bootstrap_estimates" = .check_is_bootstrap_estimate(user_input, parameter_name, additional_function_parameter),
            "delay_matrix_column_fit" = .is_value_in_accepted_values_vector(user_input, parameter_name),
            stop(paste("Checking function for type", input_type, "not found."))
    )
    }
  return(TRUE)
}
