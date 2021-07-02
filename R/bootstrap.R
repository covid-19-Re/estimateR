#' Obtain a bootstrap replicate from incidence data
#'
#' Apply a bootstrapping procedure on some original incidence data.
#' Estimating Re over many bootstrapped replicates allows one to estimate
#' the uncertainty over the estimated Re value due to observation noise.
#' In future updates, the user will be able to choose between different
#' bootstrapping procedures.
#' For now, only one bootstrapping function is implemented.
#' It performs a non-parametric block bootstrapping.
#'
#' @inheritParams module_methods
#' @inheritParams module_structure
#' @inheritDotParams .block_bootstrap -incidence_input
#'
#' @return a module output object. A boostrapped replicate. TODO improve
#' @export
get_bootstrap_replicate <- function( incidence_data,
                                     bootstrapping_method = "non-parametric block boostrap",
                                     simplify_output = TRUE,
                                     ... ) {

  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(bootstrapping_method, "bootstrapping_method"),
                                  list(simplify_output, "boolean")))



  dots_args <- .get_dots_as_list(...)
  input <- .get_module_input(incidence_data)

  if(bootstrapping_method == "non-parametric block boostrap") {
    bootstrapped_incidence <- do.call(
      '.block_bootstrap',
      c(list(incidence_input = input),
        .get_shared_args(list(.block_bootstrap,
                              .block_bootstrap_overlap_func,
                              .smooth_LOESS),
                         dots_args))
    )
  } else {
    bootstrapped_incidence <- .make_empty_module_output()
  }

  if(simplify_output) {
    bootstrapped_incidence <- .simplify_output(bootstrapped_incidence)
  }

  return(bootstrapped_incidence)
}

#' Apply block-bootstrapping procedure to module input
#'
#'\code{.block_bootstrap} returns a block-bootstrapped replicate
#'of the incidence. Incidence should be a vector of non-negative values
#'
#'This function works by resampling blocks of differences (on the log-scale)
#' between the original data and a smoothed version of the original data.
#'
#' @param round_incidence boolean. If \code{TRUE}, the bootstrapped incidence is rounded to the nearest integer.
#' @inheritParams module_methods
#' @inheritParams inner_module
#' @inheritDotParams .block_bootstrap_overlap_func -incidence_vector
#' @inheritDotParams smooth_incidence -simplify_output -incidence_data
#'
#' @return a module output object. bootstrapped incidence.
.block_bootstrap <- function(incidence_input, round_incidence = TRUE, smoothing_method = "LOESS", ...) {
  
  .are_valid_argument_values(list(list(incidence_input, "module_input"),
                                  list(smoothing_method, "smoothing_method"),
                                  list(round_incidence, "boolean")))
  

  dots_args <- .get_dots_as_list(...)

  incidence_vector <- .get_values(incidence_input)

  log_original <-log(incidence_vector + 1)

  smoothed_log <- do.call(
    'smooth_incidence',
    c(list(incidence_data = log_original,
           smoothing_method = smoothing_method),
      .get_shared_args(.smooth_LOESS, dots_args))
  )

  diff_smoothed_original <- log_original - smoothed_log

  bootstrapped_diff <- do.call(
    '.block_bootstrap_overlap_func',
    c(list(incidence_vector = diff_smoothed_original),
      .get_shared_args(.block_bootstrap_overlap_func, dots_args))
  )

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
#' @param incidence_vector numeric vector. Original incidence to bootstrap over.
#' @param block_size integer. Size of a bootstrapping block.
#' @param keep_weekdays_aligned boolean.
#' Set to \code{FALSE} if not daily incidence, or if no weekly noise pattern that would require to apply errors to the same day of the week as they were in the original data.
#'
#' @return numeric vector. Bootstrapped differences.
.block_bootstrap_overlap_func <- function(incidence_vector, block_size = 10, keep_weekdays_aligned = TRUE){

  .are_valid_argument_values(list(list(incidence_vector, "numeric_vector"),
                                  list(block_size, "positive_integer"),
                                  list(keep_weekdays_aligned, "boolean")))
  
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
