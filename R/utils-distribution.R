#TODO fill in doc
#' Title
#'
#' @param f
#' @param distribution
#'
#' @return
.get_distribution_parms <- function(distribution, f){
  # Remove the name element from the distribution list
  distribution <- within(distribution, rm(name))

  # Only keep elements of 'distribution' that are arguments of function f
  distribution_parms <- distribution[names(distribution) %in% methods::formalArgs(f)]

  return(distribution_parms)
}

#TODO fill in doc
#' Title
#'
#' @param distribution
#' @param function_prefix
#'
#' @return
.get_distribution_function <- function(distribution, function_prefix){
  f <- get(paste0(function_prefix, distribution[["name"]]), envir = loadNamespace("stats"))
  return(f)
}

#TODO fill in doc
#' Title
#'
#' @param distribution
#' @param p
#'
#' @return
.get_quantiles <- function(distribution, p) {

  q_distribution_function <- .get_distribution_function(distribution = distribution,
                                                        function_prefix = "q")

  distribution_parms <- .get_distribution_parms(distribution = distribution,
                                                f = q_distribution_function)

  return(do.call(q_distribution_function, c(list(p = p), distribution_parms)))
}

#TODO fill doc
#' Title
#'
#' @param distribution
#' @param n
#'
#' @return
.sample_from_distribution <- function(distribution, n) {
  r_distribution_function <- .get_distribution_function(distribution = distribution,
                                                        function_prefix = "r")

  distribution_parms <- .get_distribution_parms(distribution = distribution,
                                                f = r_distribution_function)

  return(do.call(r_distribution_function, c(list(n = n), distribution_parms)))
}

#TODO fill in doc
#' Title
#'
#' @param distribution
#' @param right_boundary
#' @param offset_by_one
#'
#' @return
.get_discretized_distribution <- function(distribution, right_boundary, offset_by_one){

  p_distribution_function <-  .get_distribution_function(distribution = distribution,
                                                         function_prefix = "p")

  distribution_parms <- .get_distribution_parms(f = p_distribution_function,
                                                distribution = distribution)

  if(offset_by_one) {
    x_values <- c(0, seq(from = 1.5, to = right_boundary, by = 1))
  } else {
    x_values <- c(0, seq(from = 0.5, to = right_boundary, by = 1))
  }

  cdf_values <- do.call(p_distribution_function, c(list(q=x_values), distribution_parms))

  if(length(cdf_values) == 1) {
    return(0)
  } else {
    return(diff(cdf_values))
  }
}

#' Title
#'
#' @param distribution
#' @param max_quantile
#'
#' @return
.get_right_boundary_for_distribution_vector <- function(distribution, max_quantile){
  right_boundary <- ceiling(.get_quantiles(distribution, p = max_quantile)) + 1

  # Set the right boundary to at least two
  right_boundary <- max(right_boundary, 2)

  return(right_boundary)
}

#TODO add checks for other distributions (lognormal, uniform, weibull, truncated_normal,...)
#TODO fill in doc
#' Title
#'
#' @param distribution
#'
#' @return a boolean value.
.is_valid_distribution <- function(distribution){

  if (class(distribution) != "list") {
    stop("Distribution must be a named list.")
  }

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

#TODO add details on the discretization
#TODO fill documentation
#TODO specify format of distribution: distribution <- list(name = "gamma", shape = 2, scale = 4)
#TODO test (test that vector sums up to 1)
#' Build a delay distribution vector
#'
#' @param distribution
#' @param max_quantile numeric value between 0 and 1. TODO write what max_quantile does
#' @param offset_by_one boolean. Gamma distribution comes from fit on (raw_data + 1) to accommodate zeroes in the raw_data.
#'
#' @return numeric vector.
#' @export
#'
#' @examples
#' #TODO add example
build_delay_distribution <- function(distribution,
                                     max_quantile = 0.999,
                                     offset_by_one = FALSE){

  if(!.is_valid_distribution(distribution)) {
    return(0)
  }

  right_boundary <- .get_right_boundary_for_distribution_vector(distribution = distribution,
                                                                max_quantile = max_quantile)

  distribution_vector <- .get_discretized_distribution(distribution = distribution,
                                                       right_boundary = right_boundary,
                                                       offset_by_one = offset_by_one)

  return(distribution_vector)
}
