#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom rlang :=
NULL


#' Utility functions for input validity.
#'
#' @param string_user_input A string containing the value that the user passed for the tested string type parameter.
#' @param parameter_name A string containing the name the tested parameter had in the initial function in which it was passed.
#' @param incidence_data_length A number representing the length of the given incidence data. 
#' 
#' @return TRUE if all tests were passed. Throws an error otherwise.
#' 
#' @name validation_utility_params
NULL