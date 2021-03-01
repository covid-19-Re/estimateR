
#TODO fill in documentation
#' Obtain a bootstrap replicate from incidence data
#'
#' @param original_data
#' @param bootstrapping_method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
get_bootstrap_replicate <- function( original_data, bootstrapping_method = "non-parametric block boostrap", ... ) {

  input <- get_module_input(original_data)

  bootstrapped_incidence <- dplyr::case_when(
    bootstrapping_method == "non-parametric block boostrap" ~ block_bootstrap(input, ... ),
    TRUE ~ rep(NA_real_, length.out = get_input_length(input))
  )
  # return(bootstrapped_incidence)
}



#TODO apply block-bootstrapping bug fix from Shiny dashboard code (cf commit by Jana)
#TODO test what happens with NAs in the original vector

#TODO fill in documentation
#' Apply block-bootstrapping procedure
#'
#' @param incidence_input
#' @param block_size
#' @param days_incl
#'
#' @return
#'
#' @examples
block_bootstrap <- function(incidence_input, block_size = 10, days_incl = 21) {

  incidence_vector <- incidence_input$values

  log_original <- dplyr::if_else(incidence_data !=0, log(incidence_vector + 1), 0)
  smoothed_log <- smooth_incidence(log_original, smoothing_method = "LOESS")
  diff_smoothed_original <- dplyr::if_else(log_original != 0, log_original - smoothed_log, 0)

  bootstrapped_diff <- block_bootstrap_overlap_func(diff_smoothed_original, block_size)

  bootstrapped_incidence <- exp(bootstrapped_diff + smoothed_log) -1
  bootstrapped_incidence[bootstrapped_incidence<0] <- 0
  bootstrapped_incidence <- round(bootstrapped_incidence) #TODO we don't necessarily want to round the incidence, make it an argument

  return(get_module_output(bootstrapped_incidence, incidence_input))
}

#TODO fill in function doc
#' Title
#'
#' @param incidence_vector
#' @param block_size
#' @param keep_weekdays_aligned
#'
#' @return
#' @export
#'
#' @examples
block_bootstrap_overlap_func <- function(incidence_vector, block_size = 10, keep_weekdays_aligned = TRUE){

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
