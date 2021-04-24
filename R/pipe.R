#TODO build pipe_functions that reads in "configuration" file

#TODO build way to pass extra parameters to functions inside pipe (with ...)

#TODO add a utility to summarize the uncertainty in the get_block_bootstrapped_estimate pipe

#TODO generalize get_block_bootstrapped_estimate to any bootstrapping that works the same way (parallel estimation and summary as a final step)

#TODO expand on the output
#' Estimate Re from incidence and estimate uncertainty with block-bootstrapping
#'
#' An estimation of the effective reproductive number through time is made with \code{smooth_deconvolve_estimate}
#' on the original incidence data.
#' Then, the same estimation is performed on a number of bootstrap samples built from the original incidence data.
#' The estimate on the original data is output along with confidence interval boundaries
#' built from the distribution of bootstrapped estimates.
#'
#' A base assumption made here is that the delay between infection can be split into two independent delays
#' (with the second delay being optional).
#' The first waiting time corresponds to the time between the infection event and the onset of symptoms.
#' The second (optional) waiting time corresponds to the delay between the onset of symptoms and the case observation.
#' If this second waiting time is omitted, observation is assumed to have happened at the end of the first waiting time.
#' In latter case, the first waiting time distribution does not need to correspond to an incubation period per se,
#' but it must in any case correspond to the delay between infection and case observation.
#'
#' #TODO clarify input
#'
#' @param N_bootstrap_replicates integer. Number of bootstrap samples.
#' @param uncertainty_summary_method string. 'NONE' or one of the possible strings in \code{\link{summarise_uncertainty}}
#' @inheritParams smooth_deconvolve_estimate
#'
#' @return tibble. Re estimations along with results from each pipeline step.
#' @export
get_block_bootstrapped_estimate <- function(incidence_data,
                                            N_bootstrap_replicates = 100,
                                            smoothing_method = "LOESS",
                                            deconvolution_method = "Richardson-Lucy delay distribution",
                                            estimation_method = "EpiEstim sliding window",
                                            uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                            delay_incubation,
                                            delay_onset_to_report = c(1.0),
                                            estimation_window = 3,
                                            mean_serial_interval = 4.8,
                                            std_serial_interval  = 2.3,
                                            mean_Re_prior = 1,
                                            ref_date = NULL,
                                            time_step = "day",
                                            verbose = FALSE){

  #TODO if verbose:set up

  # Display progress bar
  progress_bar <- utils::txtProgressBar(min = 0, max = N_bootstrap_replicates + 1, style = 3)
  utils::setTxtProgressBar(progress_bar, 0)

  # Prepare delay distribution vector or matrix early on as it spares the need to redo the same operation for each bootstrap replicate
  total_delay_distribution <- convolve_delay_inputs(delay_incubation,
                                                    delay_onset_to_report,
                                                    n_report_time_steps = length(incidence_data))

  original_result <- smooth_deconvolve_estimate(incidence_data,
                                                smoothing_method = smoothing_method,
                                                deconvolution_method = deconvolution_method,
                                                estimation_method = estimation_method,
                                                delay_incubation = total_delay_distribution,
                                                estimation_window = estimation_window,
                                                mean_serial_interval = mean_serial_interval,
                                                std_serial_interval  = std_serial_interval,
                                                mean_Re_prior = mean_Re_prior,
                                                ref_date = ref_date,
                                                time_step = time_step,
                                                output_Re_only = FALSE,
                                                verbose = verbose)

  original_result$bootstrap_id <- 0

  bootstrapping_results <- list(original_result)

  for(i in 1:N_bootstrap_replicates) {

    utils::setTxtProgressBar(progress_bar, i)

    bootstrapped_incidence <- get_bootstrap_replicate(incidence_data = incidence_data,
                                                      bootstrapping_method = "non-parametric block boostrap")

    bootstrapping_result <- smooth_deconvolve_estimate(bootstrapped_incidence,
                                                       smoothing_method = smoothing_method,
                                                       deconvolution_method = deconvolution_method,
                                                       estimation_method = estimation_method,
                                                       delay_incubation = total_delay_distribution,
                                                       mean_serial_interval = mean_serial_interval,
                                                       std_serial_interval  = std_serial_interval,
                                                       estimation_window = estimation_window,
                                                       ref_date = ref_date,
                                                       time_step = time_step,
                                                       output_Re_only = FALSE,
                                                       verbose = verbose)

    bootstrapping_result$bootstrap_id <- i

    bootstrapping_results <- c(bootstrapping_results, list(bootstrapping_result))
  }

  bootstrapped_estimates <- dplyr::bind_rows(bootstrapping_results)

  original_estimates <- bootstrapped_estimates %>%
    dplyr::filter(.data$bootstrap_id == 0)

  bootstrapped_estimates <- bootstrapped_estimates %>%
    dplyr::filter(.data$bootstrap_id > 0)

  estimates_with_uncertainty <- summarise_uncertainty(original_estimates = original_estimates,
                                                      bootstrapped_estimates = bootstrapped_estimates,
                                                      uncertainty_summary_method = uncertainty_summary_method,
                                                      Re_estimate_col = "R_mean",
                                                      bootstrap_id_col = "bootstrap_id",
                                                      time_step = time_step)

  # Close progress bar
  utils::setTxtProgressBar(progress_bar, N_bootstrap_replicates + 1)
  close(progress_bar)

  return(estimates_with_uncertainty)
}


#TODO rename function to remove deconvolution word
#TODO update doc
#' Estimate Re from incidence data
#'
#' This pipe function combines a smoothing step using (to remove noise from the original observations),
#' a deconvolution step (to retrieve infection events from the observed delays),
#' and an Re estimation step (wrapping around EpiEstim).
#'
#' The smoothing step uses the LOESS method by default.
#' The deconvolution step uses the Richardson-Lucy algorithm  by default.
#' The estimation uses the Cori method with a sliding window  by default.
#'
#'#TODO clarify input
#'
#' @param incidence_data numeric. TODO define acceptable format
#' @param smoothing_method string. Method used to smooth the original incidence data.
#'  Options are: "LOESS", implemented in \code{\link{smooth_LOESS}}.
#' @param deconvolution_method string. Method used to infer timings of infection
#' events from the original incidence data (aka deconvolution).
#' Options are "Richardson-Lucy delay distribution".
#' @param estimation_method string. Method used to estimate reproductive number
#' values through time from the reconstructed infection timings.
#' Options are "EpiEstim sliding window".
#' @param delay_incubation
#' @param delay_onset_to_report
#' @param estimation_window integer. Only use if \code{estimation_method = "EpiEstim sliding window"}.
#' @param mean_serial_interval numeric. see \code{\link{estimate_Re}}
#' @param std_serial_interval numeric. see \code{\link{estimate_Re}}
#' @param mean_Re_prior numeric value. Mean of prior distribution on Re.
#' @param verbose boolean. Print verbose output?
#' @param output_Re_only boolean. Should the output only contain Re estimates? (as opposed to containing results for each intermediate step)
#' @param ref_date Date. Optional. Date of the first data entry in \code{incidence_data}
#' @param time_step string. Time between two consecutive incidence datapoints. "day", "2 days", "week", "year"... (see \code{\link[base]{seq.Date}} for details)
#'
#'#TODO add details on output formatting
#' @return effective reproductive number estimates through time.
#' @export
smooth_deconvolve_estimate <- function(incidence_data,
                                       smoothing_method = "LOESS",
                                       deconvolution_method = "Richardson-Lucy delay distribution",
                                       estimation_method = "EpiEstim sliding window",
                                       delay_incubation,
                                       delay_onset_to_report = c(1.0),
                                       estimation_window = 3,
                                       mean_serial_interval = 4.8,
                                       std_serial_interval  = 2.3,
                                       mean_Re_prior = 1,
                                       output_Re_only = TRUE,
                                       ref_date = NULL,
                                       time_step = "day",
                                       verbose = FALSE) {

  smoothed_incidence <- smooth_incidence(incidence_data = incidence_data,
                                         smoothing_method = smoothing_method)

  deconvolved_incidence <- deconvolve_incidence(incidence_data = smoothed_incidence,
                                                deconvolution_method = deconvolution_method,
                                                delay_incubation = delay_incubation,
                                                delay_onset_to_report = delay_onset_to_report,
                                                start_date = ref_date,
                                                time_step = time_step,
                                                verbose = verbose)

  estimated_Re <- estimate_Re(incidence_data = deconvolved_incidence,
                              estimation_method = estimation_method,
                              estimation_window = estimation_window,
                              mean_serial_interval = mean_serial_interval,
                              std_serial_interval = std_serial_interval,
                              mean_Re_prior = mean_Re_prior)

  if(output_Re_only) {
    return(estimated_Re)
  } else {
    merged_results <- merge_outputs(
      list("observed_incidence" = incidence_data,
           "smoothed_incidence" = smoothed_incidence,
           "deconvolved_incidence" = deconvolved_incidence,
           "R_mean" = estimated_Re),
      ref_date = ref_date,
      time_step = time_step)

    return(merged_results)
  }
}
