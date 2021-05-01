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
#' @return TRUE if no error was raised.
.check_is_probability_distr_vector <- function(distribution, tolerate_NAs = FALSE, tolerance_on_sum = 1E-3) {

  .check_class(distribution, "vector", mode = "numeric")

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
#' Throws an error if not a list, or not a list with the appropriate elements.
#' Returns FALSE if parameter values return an improper distribution (if gamma distr)
#'
#' @inheritParams distribution
#'
#' @return a boolean value.
.is_valid_distribution <- function(distribution){

  .check_class(distribution, "list")

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
#' @inherit empirical_delay_data_format
#' @param delay object to be tested
#'
#' @return boolean. \code{TRUE} if the input is a dataframe in the proper format.
.check_is_empirical_delay_data <- function(delay){
  if(is.data.frame(delay)) {

    if("event_date" %!in% colnames(delay)) {
      stop("Missing 'event_date' column in dataframe.")
    }
    .check_class(delay$event_date, "Date")

    if ("report_delay" %!in% colnames(delay)) {
      stop("Missing 'report_delay' column in dataframe.")
    }
    .check_class(delay$report_delay, "vector", mode = "numeric")

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
