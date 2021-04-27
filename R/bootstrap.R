
#' Obtain a bootstrap replicate from incidence data
#'
#' Apply a bootstrapping procedure on some original incidence data.
#' Estimating Re over many bootstrapped replicates allows one to estimate
#' the uncertainty over the estimated Re value due to observation noise.
#' In future versions, the user will be able to choose between different
#' bootstrapping procedures.
#' For now, only one bootstrapping function is implemented.
#' It performs a non-parametric block bootstrapping.
#'
#' @param incidence_data TODO Format still to define.
#' @param simplify_output boolean. Return a numeric vector instead of module output object if output offset is zero.
#' @param ... Additional parameters. TODO add details
#' @inheritParams smooth_deconvolve_estimate
#'
#' @return a module output object. A boostrapped replicate.
#' @export
get_bootstrap_replicate <- function( incidence_data,
                                     bootstrapping_method = "non-parametric block boostrap",
                                     simplify_output = TRUE,
                                     ... ) {
  
  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(bootstrapping_method, "bootstrapping_method"),
                                  list(simplify_output, "boolean")))
  
  
  input <- .get_module_input(incidence_data)

  if(...length() > 0) {
    dots <- list(...)
  } else {
    dots <- list()
  }

  if(bootstrapping_method == "non-parametric block boostrap") {
    block_bootstrap_args <- names(formals(block_bootstrap))

    bootstrapped_incidence <- do.call(
      'block_bootstrap',
      c(list(incidence_input = input), dots[names(dots) %in% block_bootstrap_args])
    )
  } else {
    bootstrapped_incidence <- .make_empty_module_output()
  }

  if(simplify_output) {
    bootstrapped_incidence <- .simplify_output(bootstrapped_incidence)
  }

  return(bootstrapped_incidence)
}


#TODO polish doc (inclding details on use of LOESS)
#' Apply block-bootstrapping procedure to module input
#'
#'\code{.block_bootstrap} returns a block-bootstrapped replicate
#'of the incidence. Incidence should be a vector of non-negative values
#'
#'This function works by resampling blocks of differences (on the log-scale)
#' between the original data and a smoothed version of the original data.
#'
#' @param incidence_input module input. Original incidence to bootstrap over.
#' @param block_size integer. Size of a bootstrapping block.
#' @param round_incidence boolean. Round the bootstrapped incidence?
#' @inheritParams smooth_LOESS
#'
#' @return a module output object
#' @export
block_bootstrap <- function(incidence_input, block_size = 10, data_points_incl = 21, degree = 1, round_incidence = TRUE) {

  .are_valid_argument_values(list(list(incidence_input, "module_input"),
                                  list(block_size, "non_negative_number"),
                                  list(data_points_incl, "non_negative_number"),
                                  list(degree, "non_negative_number"),
                                  list(round_incidence, "boolean")))
  
  incidence_vector <- .get_values(incidence_input)

  log_original <-log(incidence_vector + 1)
  smoothed_log <- smooth_incidence(log_original, smoothing_method = "LOESS", data_points_incl = data_points_incl, degree = degree)
  diff_smoothed_original <- log_original - smoothed_log

  bootstrapped_diff <- .block_bootstrap_overlap_func(diff_smoothed_original, block_size)

  bootstrapped_incidence <- exp(bootstrapped_diff + smoothed_log) -1

  bootstrapped_incidence[bootstrapped_incidence<0] <- 0

  if(round_incidence) {
    bootstrapped_incidence <- round(bootstrapped_incidence)
  }

  return(.get_module_output(bootstrapped_incidence, incidence_input))
}

#' Helper function for block-bootstrapping
#'
#' Builds a bootstrapped vector of errors.
#'
#' @param incidence_vector module input. Original incidence to bootstrap over.
#' @param block_size integer. Size of a bootstrapping block.
#' @param keep_weekdays_aligned boolean. Set to FALSE if not daily incidence, or if no weekly noise pattern.
#'
#' @return numeric vector. Bootstrapped differences.
.block_bootstrap_overlap_func <- function(incidence_vector, block_size = 10, keep_weekdays_aligned = TRUE){

  bootstrapped_incidence <-c()

  if(keep_weekdays_aligned) {
    # get the weekdays for each position of incidence_vector
    weekdays_index <- (1:length(incidence_vector)) %% 7
    weekdays_index[which(weekdays_index==0)] <- 7

    last_day_index <- 7
  }

  while(length(bootstrapped_incidence) < length(incidence_vector)){
    start_index <- sample(1:(length(incidence_vector)-block_size+1), 1)
    sampled_index <- start_index:(start_index+block_size-1)

    if(keep_weekdays_aligned) {
      sampled_weekdays <- weekdays_index[sampled_index]

      # make sure the day related to the first sample is after the previous bootstrapped_incidence
      first_day_index <- which(sampled_weekdays==last_day_index)[1] + 1
      bootstrapped_incidence_index <- sampled_index[first_day_index:block_size]

      last_day_index <- utils::tail(weekdays_index[bootstrapped_incidence_index],1)
    } else {
      bootstrapped_incidence_index <- sampled_index
    }

    bootstrapped_incidence <- c(bootstrapped_incidence, incidence_vector[bootstrapped_incidence_index])
  }

  # take the same length as previous incidence_vector
  bootstrapped_incidence <- bootstrapped_incidence[1:length(incidence_vector)]

  return(bootstrapped_incidence)
}
