
#TODO fill in doc
#' Title
#'
#' @param incidence_vector
#' @param N_bootstrap_replicates
#' @param smoothing_method
#' @param shape_incubation
#' @param scale_incubation
#' @param shape_onset_to_report
#' @param scale_onset_to_report
#' @param length_out
#' @param n_random_samples
#' @param minimul_cumul_incidence
#' @param estimation_window
#' @param mean_serial_interval
#' @param std_serial_interval
#' @param mean_Re_prior
#'
#' @return
#' @export
#'
#' @examples
get_block_bootstrapped_estimate <- function(incidence_vector,
                                            N_bootstrap_replicates = 100,
                                            smoothing_method = "LOESS",
                                            shape_incubation,
                                            scale_incubation,
                                            shape_onset_to_report,
                                            scale_onset_to_report,
                                            length_out = 200,
                                            n_random_samples = 1E6,
                                            minimul_cumul_incidence = 0,
                                            estimation_window = 3,
                                            mean_serial_interval = 4.8,
                                            std_serial_interval  = 2.3,
                                            mean_Re_prior = 1){


  delay_distribution_vector <- get_vector_constant_waiting_time_distr(shape_incubation,
                                                                      scale_incubation,
                                                                      shape_onset_to_report,
                                                                      scale_onset_to_report,
                                                                      length_out = length_out,
                                                                      n_random_samples = n_random_samples)

  initial_delta <- ceiling(min(which(cumsum(delay_distribution_vector) > 0.5)) - 1)

  original_result <- smooth_deconvolve_estimate(incidence_vector,
                                                delay_distribution_vector,
                                                initial_delta,
                                                output_Re_only = FALSE) #TODO adjust input parameters

  original_result$bootstrap_id <- 0

  bootstrapping_results <- list(original_result)

  for(i in 1:N_bootstrap_replicates) {

    bootstrapped_incidence <- get_bootstrap_replicate(incidence_data = incidence_vector,
                                                      bootstrapping_method = "non-parametric block boostrap")

    bootstrapping_result <- smooth_deconvolve_estimate(bootstrapped_incidence,
                                                       delay_distribution_vector,
                                                       initial_delta,
                                                       output_Re_only = FALSE) #TODO adjust input parameters

    bootstrapping_result$bootstrap_id <- i

    bootstrapping_results <- c(bootstrapping_results, list(bootstrapping_result))
  }

  return(dplyr::bind_rows(bootstrapping_results))
}


#TODO fill in doc
#TODO rename parameters to make it more explicit where they belong
#' Title
#'
#' @param incidence_vector
#' @param delay_distribution_vector
#' @param initial_delta
#' @param minimul_cumul_incidence
#' @param estimation_window
#' @param mean_serial_interval
#' @param std_serial_interval
#' @param mean_Re_prior
#' @param time_units_in_the_past
#' @param max_iterations
#' @param verbose
#' @param output_Re_only
#' @param ref_date
#' @param time_step
#'
#' @return
#' @export
#'
#' @examples
smooth_deconvolve_estimate <- function(incidence_vector,
                                       delay_distribution_vector,
                                       initial_delta,
                                       minimul_cumul_incidence = 0,
                                       estimation_window = 3,
                                       mean_serial_interval = 4.8,
                                       std_serial_interval  = 2.3,
                                       mean_Re_prior = 1,
                                       time_units_in_the_past = 30,
                                       max_iterations = 100,
                                       verbose = FALSE,
                                       output_Re_only = TRUE,
                                       ref_date = NULL,
                                       time_step = "day") {

  smoothed_incidence <- smooth_incidence(incidence_data = incidence_vector,
                                         smoothing_method = "LOESS")

  deconvolved_incidence <- deconvolve_incidence(incidence_data = smoothed_incidence,
                                                deconvolution_method = "Richardson-Lucy delay distribution",
                                                delay_distribution_vector,
                                                initial_delta = initial_delta,
                                                time_units_in_the_past = time_units_in_the_past,
                                                max_iterations = max_iterations,
                                                verbose = verbose)

  estimated_Re <- estimate_Re(incidence_data = deconvolved_incidence,
                              estimation_method = "EpiEstim sliding window",
                              minimul_cumul_incidence = minimul_cumul_incidence,
                              estimation_window = estimation_window,
                              mean_serial_interval = mean_serial_interval,
                              std_serial_interval = std_serial_interval,
                              mean_Re_prior = mean_Re_prior)

  if(output_Re_only) {
    return(estimated_Re)
  } else {
    merged_results <- merge_outputs(
      list(observed_incidence = incidence_vector,
           smoothed_incidence = smoothed_incidence,
           deconvolved_incidence = deconvolved_incidence,
           R_mean = estimated_Re),
      ref_date = ref_date,
      time_step = time_step)

    return(merged_results)
  }

}

#TODO build pipe_functions that reads in "configuration" file
