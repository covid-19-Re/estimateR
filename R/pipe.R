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

  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(N_bootstrap_replicates, "non_negative_number"),
                                  list(smoothing_method, "smoothing_method"),
                                  list(deconvolution_method, "deconvolution_method"),
                                  list(estimation_method, "estimation_method"),
                                  list(uncertainty_summary_method, "uncertainty_summary_method"),
                                  list(delay_incubation, "delay_object", .get_input_length(incidence_data)),
                                  list(delay_onset_to_report, "delay_object", .get_input_length(incidence_data)),
                                  list(ref_date, "null_or_date"),
                                  list(time_step, "time_step")))

  dots_args <- .get_dots_as_list(...)

  index_col <- "idx"

  # Display progress bar
  # progress_bar <- utils::txtProgressBar(min = 0, max = N_bootstrap_replicates + 1, style = 3)
  # utils::setTxtProgressBar(progress_bar, 0)

  # Prepare delay distribution vector or matrix early on as it spares the need to redo the same operation for each bootstrap replicate
  total_delay_distribution <- do.call(
    'convolve_delay_inputs',
    c(list(delay_incubation = delay_incubation,
           delay_onset_to_report = delay_onset_to_report,
           n_report_time_steps = length(incidence_data),
           ref_date = ref_date,
           time_step = time_step),
      .get_shared_args(convolve_delay_inputs, dots_args))
  )

  smooth_deconvolve_estimate_dots_args <- .get_shared_args(list(.smooth_LOESS,
                                                                .deconvolve_incidence_Richardson_Lucy,
                                                                .estimate_Re_EpiEstim_sliding_window),
                                                           dots_args)

  original_result <- do.call(
    'smooth_deconvolve_estimate',
    c(list(incidence_data = incidence_data,
           smoothing_method = smoothing_method,
           deconvolution_method = deconvolution_method,
           estimation_method = estimation_method,
           delay_incubation = total_delay_distribution,
           ref_date = NULL,
           output_Re_only = FALSE,
           include_index = TRUE,
           index_col = index_col),
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
             ref_date = NULL,
             output_Re_only = FALSE,
             include_index = TRUE,
             index_col = index_col),
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
                                                      index_col = index_col)
  if(!is.null(ref_date)) {
    estimates_with_uncertainty <- .add_date_column(estimates_with_uncertainty,
                                                   ref_date = ref_date,
                                                   time_step = time_step,
                                                   index_col= index_col,
                                                   keep_index_col = FALSE)
  }

  # Close progress bar
  # utils::setTxtProgressBar(progress_bar, N_bootstrap_replicates + 1)
  # close(progress_bar)

  return(estimates_with_uncertainty)
}


#TODO rename function to remove deconvolution word
#TODO replace internals with call to get_infections_from_incidence
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
#' @inheritParams dating
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

  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(smoothing_method, "smoothing_method"),
                                  list(deconvolution_method, "deconvolution_method"),
                                  list(estimation_method, "estimation_method"),
                                  list(delay_incubation, "delay_object", .get_input_length(incidence_data)), # need to pass length of incidence data as well in order
                                  list(delay_onset_to_report, "delay_object", .get_input_length(incidence_data)), # to validate when the delay is passed as a matrix
                                  list(ref_date, "null_or_date"),
                                  list(time_step, "time_step"),
                                  list(output_Re_only, "boolean")))

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
    merged_results <- do.call(
      'merge_outputs',
      c(list(output_list = list("observed_incidence" = incidence_data,
                                "smoothed_incidence" = smoothed_incidence,
                                "deconvolved_incidence" = deconvolved_incidence,
                                "R_mean" = estimated_Re),
             ref_date = ref_date,
             time_step = time_step),
        .get_shared_args(merge_outputs, dots_args))
    )

    return(merged_results)
  }
}

#TODO doc
#TODO test
#TODO test output when output_infection_incidence_only = FALSE
#TODO rework on way delays are input (delay_incubation and delay_onset_to_report args)
#' Get timeseries of infection events from incidence data of delayed observations
#'
#' @param incidence_data
#' @param smoothing_method
#' @param deconvolution_method
#' @param delay_incubation
#' @param delay_onset_to_report
#' @param is_partially_reported_data
#' @param delay_distribution_final_report
#' @param output_infection_incidence_only
#' @inheritDotParams merge_outputs -output_list -include_index -index_col
#' @inheritDotParams correct_for_partially_observed_data -incidence_data -delay_distribution_final_report
#'
#' @return
#' @export
get_infections_from_incidence <- function(incidence_data,
                                          smoothing_method = "LOESS",
                                          deconvolution_method = "Richardson-Lucy delay distribution",
                                          delay_incubation,
                                          delay_onset_to_report = c(1.0),
                                          is_partially_reported_data = FALSE,
                                          delay_distribution_final_report = NULL,
                                          output_infection_incidence_only = TRUE,
                                          ...){

  .are_valid_argument_values(list(list(incidence_data, "module_input"),
                                  list(smoothing_method, "smoothing_method"),
                                  list(deconvolution_method, "deconvolution_method"),
                                  list(delay_incubation, "delay_object", .get_input_length(incidence_data)), # need to pass length of incidence data as well in order
                                  list(delay_onset_to_report, "delay_object", .get_input_length(incidence_data)), # to validate when the delay is passed as a matrix
                                  list(is_partially_reported_data, "boolean"),
                                  list(output_infection_incidence_only, "boolean")))

  dots_args <- .get_dots_as_list(...)

  original_incidence_data <- incidence_data

  if(is_partially_reported_data) {
    .are_valid_argument_values(list(list(delay_distribution_final_report, "delay_object", .get_input_length(incidence_data))))

    incidence_data <- do.call(
      'correct_for_partially_observed_data',
      c(list(incidence_data = incidence_data,
             delay_distribution_final_report = delay_distribution_final_report),
        .get_shared_args(correct_for_partially_observed_data, dots_args))
    )
  }

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

  if(output_infection_incidence_only) {
    return(deconvolved_incidence)
  } else {
    if(is_partially_reported_data) {
      output_list <- list("observed_incidence" = original_incidence_data,
                          "corrected_incidence" = incidence_data,
                          "smoothed_incidence" = smoothed_incidence,
                          "deconvolved_incidence" = deconvolved_incidence)
    } else {
      output_list <- list("observed_incidence" = incidence_data,
                          "smoothed_incidence" = smoothed_incidence,
                          "deconvolved_incidence" = deconvolved_incidence)
    }

    merged_results <- do.call(
      'merge_outputs',
      c(list(output_list = output_list),
        .get_shared_args(merge_outputs, dots_args))
    )

    return(merged_results)
  }
}



#TODO doc
#' Title
#'
#' @param partially_delayed_incidence
#' @param fully_delayed_incidence
#' @param smoothing_method
#' @param deconvolution_method
#' @param estimation_method
#' @param delay_until_partial
#' @param delay_from_partial_to_full
#' @param partial_observation_requires_full_observation
#' @param ref_date
#' @param time_step
#' @param output_Re_only
#' @param ...
#'
#' @return
#' @export
estimate_from_combined_observations <- function(partially_delayed_incidence,
                                                fully_delayed_incidence,
                                                smoothing_method = "LOESS",
                                                deconvolution_method = "Richardson-Lucy delay distribution",
                                                estimation_method = "EpiEstim sliding window",
                                                delay_until_partial,
                                                delay_from_partial_to_full,
                                                partial_observation_requires_full_observation = TRUE,
                                                ref_date = NULL,
                                                time_step = "day",
                                                output_Re_only = TRUE,
                                                ...) {

  .are_valid_argument_values(list(list(partially_delayed_incidence, "module_input"),
                                  list(fully_delayed_incidence, "module_input"),
                                  list(smoothing_method, "smoothing_method"),
                                  list(deconvolution_method, "deconvolution_method"),
                                  list(estimation_method, "estimation_method"),
                                  #TODO figure out what we do for making sure the two traces are the same size in input.
                                  list(delay_until_partial, "delay_object", .get_input_length(delay_until_partial)), # need to pass length of incidence data as well in order
                                  list(delay_from_partial_to_full, "delay_object", .get_input_length(delay_from_partial_to_full)), # to validate when the delay is passed as a matrix
                                  list(partial_observation_requires_full_observation, "boolean"),
                                  list(ref_date, "null_or_date"),
                                  list(time_step, "time_step"),
                                  list(output_Re_only, "boolean")))

  #TODO add '...' args
  infections_from_partially_delayed_observations <- get_infections_from_incidence(partially_delayed_incidence,
                                                                                  smoothing_method = smoothing_method,
                                                                                  deconvolution_method = deconvolution_method,
                                                                                  delay_incubation = delay_until_partial,
                                                                                  is_partially_reported_data = partial_observation_requires_full_observation,
                                                                                  delay_distribution_final_report = delay_from_partial_to_full,
                                                                                  output_infection_incidence_only = TRUE)
  #TODO add '...' args
  infections_from_fully_delayed_observations <- get_infections_from_incidence(fully_delayed_incidence,
                                                                              smoothing_method = smoothing_method,
                                                                              deconvolution_method = deconvolution_method,
                                                                              delay_incubation = delay_until_partial,
                                                                              delay_onset_to_report = delay_from_partial_to_full,
                                                                              is_partially_reported_data = FALSE,
                                                                              output_infection_incidence_only = TRUE)

  #TODO code inner_addition
  all_infection_events <- inner_addition(infections_from_partially_delayed_observations, infections_from_fully_delayed_observations)

  estimated_Re <- do.call(
    'estimate_Re',
    c(list(incidence_data = all_infection_events,
           estimation_method = estimation_method),
      .get_shared_args(.estimate_Re_EpiEstim_sliding_window, dots_args))
  )

  if(output_Re_only) {
    return(estimated_Re)
  } else {
    merged_results <- do.call(
      'merge_outputs',
      c(list(output_list = list("partially_delayed_observations" = partially_delayed_incidence,
                                "fully_delayed_observations" = fully_delayed_incidence,
                                "deconvolved_incidence" = all_infection_events,
                                "R_mean" = estimated_Re),
             ref_date = ref_date,
             time_step = time_step),
        .get_shared_args(merge_outputs, dots_args))
    )

    return(merged_results)
  }
}

#TODO write build_bootstrap_estimates pipe.
#TODO write pipe for combining bootstrapping estimation from combined observations.




