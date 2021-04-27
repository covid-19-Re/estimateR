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
#'
#' @param N_bootstrap_replicates integer. Number of bootstrap samples.
#' @inheritParams pipe_params
#' @inheritParams module_methods
#' @inheritParams module_structure
#' @inheritParams universal_params
#' @inheritParams delay_high
#' @inheritParams dating
#' @inheritDotParams .smooth_LOESS -incidence_input
#' @inheritDotParams .deconvolve_incidence_Richardson_Lucy -incidence_input
#' @inheritDotParams .estimate_Re_EpiEstim_sliding_window -incidence_input
#'
#' @return tibble. Re estimations along with results from each pipeline step. TODO add details
#' @export
get_block_bootstrapped_estimate <- function(incidence_data,
                                            N_bootstrap_replicates = 100,
                                            smoothing_method = "LOESS",
                                            deconvolution_method = "Richardson-Lucy delay distribution",
                                            estimation_method = "EpiEstim sliding window",
                                            uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                            delay_incubation,
                                            delay_onset_to_report = c(1.0),
                                            ref_date = NULL,
                                            time_step = "day",
                                            ...){

  dots_args <- .get_dots_as_list(...)

  # Display progress bar
  # progress_bar <- utils::txtProgressBar(min = 0, max = N_bootstrap_replicates + 1, style = 3)
  # utils::setTxtProgressBar(progress_bar, 0)

  # Prepare delay distribution vector or matrix early on as it spares the need to redo the same operation for each bootstrap replicate
  total_delay_distribution <- do.call(
    'convolve_delay_inputs',
    c(list(delay_incubation = delay_incubation,
           delay_onset_to_report = delay_onset_to_report,
           n_report_time_steps = length(incidence_data),
           start_date = ref_date,
           time_step = time_step),
      .get_shared_args(convolve_delay_inputs, dots_args))
  )

  smooth_deconvolve_estimate_dots_args <- .get_shared_args(list(.smooth_LOESS,
                                                                .deconvolve_incidence_Richardson_Lucy,
                                                                .estimate_Re_EpiEstim_sliding_window),
                                                           dots_args)

  #TODO continue here
  original_result <- do.call(
    'smooth_deconvolve_estimate',
    c(list(incidence_data = incidence_data,
           smoothing_method = smoothing_method,
           deconvolution_method = deconvolution_method,
           estimation_method = estimation_method,
           delay_incubation = total_delay_distribution,
           ref_date = ref_date,
           time_step = time_step,
           output_Re_only = FALSE),
      smooth_deconvolve_estimate_dots_args)
  )

  original_result$bootstrap_id <- 0

  bootstrapping_results <- list(original_result)

  for(i in 1:N_bootstrap_replicates) {

    # utils::setTxtProgressBar(progress_bar, i)

    bootstrapped_incidence <- do.call(
      'get_bootstrap_replicate',
      c(list(incidence_data = incidence_data,
             bootstrapping_method = "non-parametric block boostrap"),
        .get_shared_args(list(.block_bootstrap,
                              .block_bootstrap_overlap_func,
                              .smooth_LOESS),
                         dots_args))
    )

    bootstrapping_result <- do.call(
      'smooth_deconvolve_estimate',
      c(list(incidence_data = bootstrapped_incidence,
             smoothing_method = smoothing_method,
             deconvolution_method = deconvolution_method,
             estimation_method = estimation_method,
             delay_incubation = total_delay_distribution,
             ref_date = ref_date,
             time_step = time_step,
             output_Re_only = FALSE),
        smooth_deconvolve_estimate_dots_args)
    )

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
  # utils::setTxtProgressBar(progress_bar, N_bootstrap_replicates + 1)
  # close(progress_bar)

  return(estimates_with_uncertainty)
}


#TODO rename function to remove deconvolution word
#' Estimate Re from incidence data
#'
#' This pipe function combines a smoothing step using (to remove noise from the original observations),
#' a deconvolution step (to retrieve infection events from the observed delays),
#' and an Re estimation step wrapping around \code{\link[EpiEstim]{estimate_R}}.
#'
#' The \code{\link[=smooth_incidence]{smoothing step}} uses the LOESS method by default.
#' The \code{\link[=deconvolve_incidence]{deconvolution step}} uses the Richardson-Lucy algorithm  by default.
#' The \code{\link[=estimate_Re]{Re estimation}} uses the Cori method with a sliding window  by default.
#'
#' @inheritParams module_structure
#' @inheritParams module_methods
#' @inheritParams universal_params
#' @inheritParams pipe_params
#' @inheritParams delay_high
#' @inheritDotParams .smooth_LOESS -incidence_input
#' @inheritDotParams .deconvolve_incidence_Richardson_Lucy -incidence_input
#' @inheritDotParams .estimate_Re_EpiEstim_sliding_window -incidence_input
#'
#' @return effective reproductive number estimates through time TODO add details
#' @export
smooth_deconvolve_estimate <- function(incidence_data,
                                       smoothing_method = "LOESS",
                                       deconvolution_method = "Richardson-Lucy delay distribution",
                                       estimation_method = "EpiEstim sliding window",
                                       delay_incubation,
                                       delay_onset_to_report = c(1.0),
                                       ref_date = NULL,
                                       time_step = "day",
                                       output_Re_only = TRUE,
                                       ...) {

  dots_args <- .get_dots_as_list(...)

  smoothed_incidence <- do.call(
    'smooth_incidence',
    c(list(incidence_data = incidence_data,
           smoothing_method = smoothing_method),
      .get_shared_args(.smooth_LOESS, dots_args))
  )

  #TODO when generalizing to list of delays,
  # possibly take the convolution out of the deconvolution step.
  # to improve readability
  deconvolved_incidence <- do.call(
    'deconvolve_incidence',
    c(list(incidence_data = smoothed_incidence,
           deconvolution_method = deconvolution_method,
           delay_incubation = delay_incubation,
           delay_onset_to_report = delay_onset_to_report),
      .get_shared_args(list(.deconvolve_incidence_Richardson_Lucy,
                            convolve_delay_inputs),
                       dots_args))
  )

  estimated_Re <- do.call(
    'estimate_Re',
    c(list(incidence_data = deconvolved_incidence,
           estimation_method = estimation_method),
      .get_shared_args(.estimate_Re_EpiEstim_sliding_window, dots_args))
  )

  if(output_Re_only) {
    return(estimated_Re)
  } else {
    merging_args <- names(formals(merge_outputs))

    merged_results <- merge_outputs(
      output_list = list("observed_incidence" = incidence_data,
                         "smoothed_incidence" = smoothed_incidence,
                         "deconvolved_incidence" = deconvolved_incidence,
                         "R_mean" = estimated_Re),
      ref_date = ref_date,
      time_step = time_step)

    return(merged_results)
  }
}
