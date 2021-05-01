#' Get relevant parameters from distribution
#'
#' @param f distribution function from stats package.
#' @inheritParams distribution
#'
#' @return list containing elements of \code{distribution} that overlap with arguments of \code{f}
.get_distribution_parms <- function(distribution, f){
  # Remove the name element from the distribution list
  distribution <- within(distribution, rm("name"))

  # Only keep elements of 'distribution' that are arguments of function f
  distribution_parms <- distribution[names(distribution) %in% methods::formalArgs(f)]

  return(distribution_parms)
}

#' Get distribution function
#'
#' @param function_prefix character. 'd', 'q', 'p' or 'r'. see \code{\link[stats:Distributions]{stats::Distributions}}
#' @inheritParams distribution
#'
#' @return Density, distribution function, quantile function or random generation function for \code{distribution}
.get_distribution_function <- function(distribution, function_prefix){
  f <- get(paste0(function_prefix, distribution[["name"]]), envir = loadNamespace("stats"))
  return(f)
}

#' Obtain quantile values for a distribution
#'
#' @inheritParams distribution
#' @param p vector of probabilities
#'
#' @return vector of quantiles
.get_quantiles <- function(distribution, p) {

  q_distribution_function <- .get_distribution_function(distribution = distribution,
                                                        function_prefix = "q")

  distribution_parms <- .get_distribution_parms(distribution = distribution,
                                                f = q_distribution_function)

  return(do.call(q_distribution_function, c(list(p = p), distribution_parms)))
}

#' Draw samples from a probability distribution.
#'
#' @inheritParams distribution
#' @param n integer. Number of samples to draw.
#'
#' @return vector containing \code{n} samples of \code{distribution}
.sample_from_distribution <- function(distribution, n) {
  r_distribution_function <- .get_distribution_function(distribution = distribution,
                                                        function_prefix = "r")

  distribution_parms <- .get_distribution_parms(distribution = distribution,
                                                f = r_distribution_function)

  return(do.call(r_distribution_function, c(list(n = n), distribution_parms)))
}

#' Discretize a probability distribution.
#'
#' @param right_boundary positive numeric value.
#' Maximum number of time steps to discretize the \code{distribution} over.
#' @param offset_by_one boolean.
#' Set to TRUE if \code{distribution} represent the fit of data that was offset by one
#' (\code{fitted_data = original_data + 1}) to accommodate zeroes in \code{original_data}.
#' @inheritParams distribution
#'
#' @return vector containing weights of the discretized probability distribution.
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

#' Get the number of steps before a quantile is reached.
#'
#' @inheritParams distribution
#' @param max_quantile numeric value. Upper quantile that needs to be reached.
#'
#' @return number of time steps required to for the probability distribution to reach \code{max_quantile}
.get_right_boundary_for_distribution_vector <- function(distribution, max_quantile){
  right_boundary <- ceiling(.get_quantiles(distribution, p = max_quantile)) + 1

  # Set the right boundary to at least two
  right_boundary <- max(right_boundary, 2)

  return(right_boundary)
}

#TODO add details on the discretization
#' Build a discretized probability distribution vector from a delay distribution
#'
#' @inheritParams distribution
#' @inheritParams .get_discretized_distribution
#' @param max_quantile numeric value between 0 and 1. TODO add details
#' Upper quantile reached by the last element in the discretized distribution vector.
#'
#' @return vector containing the discretized probability distribution vector of \code{distribution}.
#'
#' @export
build_delay_distribution <- function(distribution,
                                     max_quantile = 0.999,
                                     offset_by_one = FALSE){

  #TODO revisit this return of false, maybe if fails we throw error always
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

#' Return probability distribution vector
#'
#' Can take a \code{distribution} list or a probability distribution vector.
#' In the first case, this function builds and return the vector of discretized probability distribution.
#' In the second case, it checks that \code{delay} is a valid discretized probability distribution and returns it.
#' See \code{\link{build_delay_distribution}} for details on the \code{distribution} list format.
#'
#' @param delay list or vector. Delay distribution to transform or validate
#' into a vector of discretized probability distribution.
#' @inheritDotParams build_delay_distribution -distribution
#'
#' @return vector of discretized probability distribution.
.get_delay_distribution <- function(delay,
                                    ...){
  if( .is_numeric_vector(delay) ) {
    .check_is_probability_distr_vector(delay)
    return(delay)
  } else if( is.list(delay) ) {
    return(build_delay_distribution(delay, ...))
  } else {
    stop("Input must either be a vector representing the discretized delay distribution or a distribution object.")
  }
}



