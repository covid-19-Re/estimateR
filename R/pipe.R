# TODO reduce duplication between bootstrapping pipe functions

# TODO polish doc (e.g. '...' args)
# TODO redoc (combine_bootstrap_and_estimation_uncertainties)
#' Estimate Re from incidence and estimate uncertainty with block-bootstrapping
#'
#' An estimation of the effective reproductive number through time is made
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
#' @inheritParams pipe_params
#' @inheritParams bootstrap_params
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
                                            combine_bootstrap_and_estimation_uncertainties = FALSE,
                                            delay_incubation,
                                            delay_onset_to_report = c(1.0),
                                            ref_date = NULL,
                                            time_step = "day",
                                            output_Re_only = TRUE,
                                            ...) {
  .are_valid_argument_values(list(
    list(incidence_data, "module_input"),
    list(N_bootstrap_replicates, "non_negative_number"),
    list(smoothing_method, "smoothing_method"),
    list(deconvolution_method, "deconvolution_method"),
    list(estimation_method, "estimation_method"),
    list(uncertainty_summary_method, "uncertainty_summary_method"),
    list(combine_bootstrap_and_estimation_uncertainties, "boolean"),
    list(delay_incubation, "delay_object", .get_input_length(incidence_data)),
    list(delay_onset_to_report, "delay_object", .get_input_length(incidence_data)),
    list(ref_date, "null_or_date"),
    list(time_step, "time_step"),
    list(output_Re_only, "boolean")
  ))

  dots_args <- .get_dots_as_list(...)

  index_col <- "idx"
  bootstrap_id_col <- "bootstrap_id"

  # Display progress bar
  # progress_bar <- utils::txtProgressBar(min = 0, max = N_bootstrap_replicates + 1, style = 3)
  # utils::setTxtProgressBar(progress_bar, 0)

  # Prepare delay distribution vector or matrix early on as it spares the need to redo the same operation for each bootstrap replicate
  total_delay_distribution <- do.call(
    "convolve_delay_inputs",
    c(
      list(
        delay_incubation = delay_incubation,
        delay_onset_to_report = delay_onset_to_report,
        n_report_time_steps = length(incidence_data),
        ref_date = ref_date,
        time_step = time_step
      ),
      .get_shared_args(convolve_delay_inputs, dots_args)
    )
  )

  smooth_deconvolve_estimate_dots_args <- .get_shared_args(
    list(
      .smooth_LOESS,
      .deconvolve_incidence_Richardson_Lucy,
      .estimate_Re_EpiEstim_sliding_window,
      .estimate_Re_EpiEstim_piecewise_constant,
      get_matrix_from_empirical_delay_distr,
      build_delay_distribution
    ),
    dots_args
  )

  original_result <- do.call(
    "estimate_Re_from_noisy_delayed_incidence",
    c(
      list(
        incidence_data = incidence_data,
        smoothing_method = smoothing_method,
        deconvolution_method = deconvolution_method,
        estimation_method = estimation_method,
        delay_incubation = total_delay_distribution,
        ref_date = NULL,
        output_Re_only = FALSE,
        include_index = TRUE,
        index_col = index_col,
        output_HPD = combine_bootstrap_and_estimation_uncertainties
      ),
      smooth_deconvolve_estimate_dots_args
    )
  )

  if (combine_bootstrap_and_estimation_uncertainties) {
    # Keep the Re estimation HPDs for later
    Re_HPDs <- original_result %>%
      dplyr::select(.data[[index_col]], .data$Re_highHPD, .data$Re_lowHPD)

    # Remove them from the original_result variable
    original_result <- original_result %>%
      dplyr::select(!c(.data$Re_highHPD, .data$Re_lowHPD))
  } else {
    Re_HPDs <- NULL
  }

  original_result[[bootstrap_id_col]] <- 0

  bootstrapping_results <- list(original_result)

  for (i in 1:N_bootstrap_replicates) {

    # utils::setTxtProgressBar(progress_bar, i)

    bootstrapped_incidence <- do.call(
      "get_bootstrap_replicate",
      c(
        list(
          incidence_data = incidence_data,
          bootstrapping_method = "non-parametric block boostrap"
        ),
        .get_shared_args(
          list(
            .block_bootstrap,
            .block_bootstrap_overlap_func,
            .smooth_LOESS
          ),
          dots_args
        )
      )
    )

    bootstrapping_result <- do.call(
      "estimate_Re_from_noisy_delayed_incidence",
      c(
        list(
          incidence_data = bootstrapped_incidence,
          smoothing_method = smoothing_method,
          deconvolution_method = deconvolution_method,
          estimation_method = estimation_method,
          delay_incubation = total_delay_distribution,
          ref_date = NULL,
          output_Re_only = FALSE,
          include_index = TRUE,
          index_col = index_col,
          output_HPD = FALSE
        ),
        smooth_deconvolve_estimate_dots_args
      )
    )

    bootstrapping_result[[bootstrap_id_col]] <- i

    bootstrapping_results <- c(bootstrapping_results, list(bootstrapping_result))
  }

  bootstrapped_estimates <- dplyr::bind_rows(bootstrapping_results)

  original_estimates <- bootstrapped_estimates %>%
    dplyr::filter(.data[[bootstrap_id_col]] == 0)

  bootstrapped_estimates <- bootstrapped_estimates %>%
    dplyr::filter(.data[[bootstrap_id_col]] > 0)

  estimates_with_uncertainty <- do.call(
    "do_uncertainty_summary",
    c(
      list(
        original_values = original_estimates,
        bootstrapped_values = bootstrapped_estimates,
        uncertainty_summary_method = uncertainty_summary_method,
        value_col = "Re_estimate",
        bootstrap_id_col = bootstrap_id_col,
        index_col = index_col,
        output_Re_only = output_Re_only,
        combine_bootstrap_and_estimation_uncertainties = combine_bootstrap_and_estimation_uncertainties,
        Re_HPDs = Re_HPDs
      ),
      .get_shared_args(.summarise_CI_bootstrap, dots_args)
    )
  )


  if (!is.null(ref_date)) {
    estimates_with_uncertainty <- .add_date_column(estimates_with_uncertainty,
      ref_date = ref_date,
      time_step = time_step,
      index_col = index_col,
      keep_index_col = FALSE
    )
  }

  # Close progress bar
  # utils::setTxtProgressBar(progress_bar, N_bootstrap_replicates + 1)
  # close(progress_bar)

  pretty_results <- do.call(
    ".prettify_result",
    c(list(data = estimates_with_uncertainty),
      .get_shared_args(.prettify_result, dots_args))
  )

  return(pretty_results)
}


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
estimate_Re_from_noisy_delayed_incidence <- function(incidence_data,
                                                     smoothing_method = "LOESS",
                                                     deconvolution_method = "Richardson-Lucy delay distribution",
                                                     estimation_method = "EpiEstim sliding window",
                                                     delay_incubation,
                                                     delay_onset_to_report = c(1.0),
                                                     ref_date = NULL,
                                                     time_step = "day",
                                                     output_Re_only = TRUE,
                                                     ...) {
  .are_valid_argument_values(list(
    list(incidence_data, "module_input"),
    list(smoothing_method, "smoothing_method"),
    list(deconvolution_method, "deconvolution_method"),
    list(estimation_method, "estimation_method"),
    list(delay_incubation, "delay_object", .get_input_length(incidence_data)), # need to pass length of incidence data as well in order
    list(delay_onset_to_report, "delay_object", .get_input_length(incidence_data)), # to validate when the delay is passed as a matrix
    list(ref_date, "null_or_date"),
    list(time_step, "time_step"),
    list(output_Re_only, "boolean")
  ))

  dots_args <- .get_dots_as_list(...)

  get_infections_from_incidence

  smoothed_incidence <- do.call(
    "smooth_incidence",
    c(
      list(
        incidence_data = incidence_data,
        smoothing_method = smoothing_method
      ),
      .get_shared_args(.smooth_LOESS, dots_args)
    )
  )

  # TODO when generalizing to list of delays,
  # possibly take the convolution out of the deconvolution step.
  # to improve readability
  deconvolved_incidence <- do.call(
    "deconvolve_incidence",
    c(
      list(
        incidence_data = smoothed_incidence,
        deconvolution_method = deconvolution_method,
        delay_incubation = delay_incubation,
        delay_onset_to_report = delay_onset_to_report
      ),
      .get_shared_args(
        list(
          .deconvolve_incidence_Richardson_Lucy,
          convolve_delay_inputs
        ),
        dots_args
      )
    )
  )

  estimated_Re <- do.call(
    "estimate_Re",
    c(
      list(
        incidence_data = deconvolved_incidence,
        estimation_method = estimation_method,
        simplify_output = FALSE
      ),
      .get_shared_args(.estimate_Re_EpiEstim_sliding_window, dots_args)
    )
  )

  if (output_Re_only) {
    merged_results <- estimated_Re
  } else {
    # TODO simplify this call (create util function)
    test_if_single_output <- try(.is_valid_module_input(estimated_Re, "estimated_Re"), silent = TRUE)
    if (!("try-error" %in% class(test_if_single_output))) {
      estimated_Re <- list("Re_estimate" = estimated_Re)
    }

    merged_results <- do.call(
      "merge_outputs",
      c(
        list(
          output_list = c(
            list(
              "observed_incidence" = incidence_data,
              "smoothed_incidence" = smoothed_incidence,
              "deconvolved_incidence" = deconvolved_incidence
            ),
            estimated_Re
          ),
          ref_date = ref_date,
          time_step = time_step
        ),
        .get_shared_args(merge_outputs, dots_args)
      )
    )
  }

  pretty_results <- do.call(
    ".prettify_result",
    c(list(data = merged_results),
      .get_shared_args(.prettify_result, dots_args))
  )

  return(pretty_results)
}


# TODO polish doc
# TODO test
#' Infer timeseries of infection events from incidence data of delayed observations
#'
#' This function takes as input incidence data of delayed observations of infections events,
#' as well as the probability distribution(s) of the delay(s).
#' It returns an inferred incidence of infection events.
#' This function can account for the observations being dependent on future delayed observations.
#' For instance, if the incidence data represents symptom onset events, usually these events
#' are dependent on a secondary delayed observation: a case confirmation typically, or
#' a hospital admission or any other type of event,
#'
#' @inheritParams module_structure
#' @inheritParams module_methods
#' @inheritParams pipe_params
#' @inheritParams delay_high
#' @inheritParams dating
#' @inheritDotParams .smooth_LOESS -incidence_input
#' @inheritDotParams .deconvolve_incidence_Richardson_Lucy -incidence_input
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
                                          ref_date = NULL,
                                          time_step = "day",
                                          ...) {
  .are_valid_argument_values(list(
    list(incidence_data, "module_input"),
    list(smoothing_method, "smoothing_method"),
    list(deconvolution_method, "deconvolution_method"),
    list(delay_incubation, "delay_object", .get_input_length(incidence_data)), # need to pass length of incidence data as well in order
    list(delay_onset_to_report, "delay_object", .get_input_length(incidence_data)), # to validate when the delay is passed as a matrix
    list(is_partially_reported_data, "boolean"),
    list(output_infection_incidence_only, "boolean"),
    list(ref_date, "null_or_date"),
    list(time_step, "time_step")
  ))

  dots_args <- .get_dots_as_list(...)

  original_incidence_data <- incidence_data

  if (is_partially_reported_data) {
    .are_valid_argument_values(list(list(delay_distribution_final_report, "delay_object", .get_input_length(incidence_data))))

    incidence_data <- do.call(
      "correct_for_partially_observed_data",
      c(
        list(
          incidence_data = incidence_data,
          delay_distribution_final_report = delay_distribution_final_report,
          ref_date = NULL,
          time_step = "day"
        ),
        .get_shared_args(correct_for_partially_observed_data, dots_args)
      )
    )
  }

  smoothed_incidence <- do.call(
    "smooth_incidence",
    c(
      list(
        incidence_data = incidence_data,
        smoothing_method = smoothing_method
      ),
      .get_shared_args(.smooth_LOESS, dots_args)
    )
  )

  # TODO when generalizing to list of delays,
  # possibly take the convolution out of the deconvolution step.
  # to improve readability
  deconvolved_incidence <- do.call(
    "deconvolve_incidence",
    c(
      list(
        incidence_data = smoothed_incidence,
        deconvolution_method = deconvolution_method,
        delay_incubation = delay_incubation,
        delay_onset_to_report = delay_onset_to_report
      ),
      .get_shared_args(
        list(
          .deconvolve_incidence_Richardson_Lucy,
          convolve_delay_inputs
        ),
        dots_args
      )
    )
  )

  if (output_infection_incidence_only) {
    merged_results <- deconvolved_incidence
  } else {
    if (is_partially_reported_data) {
      output_list <- list(
        "observed_incidence" = original_incidence_data,
        "corrected_incidence" = incidence_data,
        "smoothed_incidence" = smoothed_incidence,
        "deconvolved_incidence" = deconvolved_incidence
      )
    } else {
      output_list <- list(
        "observed_incidence" = incidence_data,
        "smoothed_incidence" = smoothed_incidence,
        "deconvolved_incidence" = deconvolved_incidence
      )
    }

    merged_results <- do.call(
      "merge_outputs",
      c(
        list(
          output_list = output_list,
          ref_date = ref_date,
          time_step = time_step
        ),
        .get_shared_args(merge_outputs, dots_args)
      )
    )
  }

  pretty_results <- do.call(
    ".prettify_result",
    c(list(data = merged_results),
      .get_shared_args(.prettify_result, dots_args))
  )

  return(pretty_results)
}


# TODO refactor arguments: we are using delay_from_partial_to_full and delay_distribution_final_report to refer to the same thing
# TODO allow to pass a third argument for delays: the convolution of all delays (to speed up bootstrapping)
# TODO polish doc
# TODO test
#' Estimate Re from delayed observations of infection events.
#'
#' This function allows for combining two different incidence timeseries.
#' The two timeseries can represent events that are differently delayed from the original infection events.
#' The two data sources must not have any overlap in the events recorded.
#' The function can account for the one of the two types of events to require
#' the future observation of the other type of event.
#' For instance, one type can be events of symptom onset, and the other be case confirmation.
#' Typically, the recording of a symptom onset event will require a future case confirmation.
#' If so, the \code{partial_observation_requires_full_observation} flag should be set to \code{TRUE}.
#'
#'
#' @inheritParams module_structure
#' @inheritParams module_methods
#' @inheritParams pipe_params
#' @inheritParams delay_high
#' @inheritParams dating
#' @inheritDotParams .smooth_LOESS -incidence_input
#' @inheritDotParams .deconvolve_incidence_Richardson_Lucy -incidence_input
#' @inheritDotParams .estimate_Re_EpiEstim_sliding_window -incidence_input
#' @inheritDotParams merge_outputs -output_list -include_index -index_col
#' @inheritDotParams correct_for_partially_observed_data -incidence_data -delay_distribution_final_report
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
  .are_valid_argument_values(list(
    list(partially_delayed_incidence, "module_input"),
    list(fully_delayed_incidence, "module_input"),
    list(smoothing_method, "smoothing_method"),
    list(deconvolution_method, "deconvolution_method"),
    list(estimation_method, "estimation_method"),
    # TODO figure out what we do for making sure the two traces are the same size in input.
    list(delay_until_partial, "delay_object", .get_input_length(partially_delayed_incidence)), # need to pass length of incidence data as well in order
    list(delay_from_partial_to_full, "delay_object", .get_input_length(fully_delayed_incidence)), # to validate when the delay is passed as a matrix
    list(partial_observation_requires_full_observation, "boolean"),
    list(ref_date, "null_or_date"),
    list(time_step, "time_step"),
    list(output_Re_only, "boolean")
  ))

  dots_args <- .get_dots_as_list(...)

  infections_from_partially_delayed_observations <- do.call(
    "get_infections_from_incidence",
    c(
      list(partially_delayed_incidence,
        smoothing_method = smoothing_method,
        deconvolution_method = deconvolution_method,
        delay_incubation = delay_until_partial,
        is_partially_reported_data = partial_observation_requires_full_observation,
        delay_distribution_final_report = delay_from_partial_to_full,
        output_infection_incidence_only = TRUE
      ),
      .get_shared_args(
        list(
          .deconvolve_incidence_Richardson_Lucy,
          convolve_delay_inputs
        ),
        dots_args
      )
    )
  )

  infections_from_fully_delayed_observations <- do.call(
    "get_infections_from_incidence",
    c(
      list(fully_delayed_incidence,
        smoothing_method = smoothing_method,
        deconvolution_method = deconvolution_method,
        delay_incubation = delay_until_partial,
        delay_onset_to_report = delay_from_partial_to_full,
        is_partially_reported_data = FALSE,
        output_infection_incidence_only = TRUE
      ),
      .get_shared_args(
        list(
          .deconvolve_incidence_Richardson_Lucy,
          convolve_delay_inputs
        ),
        dots_args
      )
    )
  )

  all_infection_events <- inner_addition(infections_from_partially_delayed_observations, infections_from_fully_delayed_observations)

  estimated_Re <- do.call(
    "estimate_Re",
    c(
      list(
        incidence_data = all_infection_events,
        estimation_method = estimation_method
      ),
      .get_shared_args(.estimate_Re_EpiEstim_sliding_window, dots_args)
    )
  )

  if (output_Re_only) {
    merged_results <- estimated_Re
  } else {
    test_if_single_output <- try(.is_valid_module_input(estimated_Re, "estimated_Re"), silent = TRUE)
    if (!("try-error" %in% class(test_if_single_output))) {
      estimated_Re <- list("Re_estimate" = estimated_Re)
    }

    merged_results <- do.call(
      "merge_outputs",
      c(
        list(
          output_list = c(
            list(
              "partially_delayed_observations" = partially_delayed_incidence,
              "fully_delayed_observations" = fully_delayed_incidence,
              "combined_deconvolved_incidence" = all_infection_events
            ),
            estimated_Re
          ),
          ref_date = ref_date,
          time_step = time_step
        ),
        .get_shared_args(merge_outputs, dots_args)
      )
    )
  }
  pretty_results <- do.call(
    ".prettify_result",
    c(list(data = merged_results),
      .get_shared_args(.prettify_result, dots_args))
  )

  return(pretty_results)
}

# TODO polish doc (e.g. '...' args)
# TODO doc
# TODO continue here
#' Estimate Re from incidence and estimate uncertainty by bootstrapping
#'
#' @inheritParams module_structure
#' @inheritParams module_methods
#' @inheritParams pipe_params
#' @inheritParams bootstrap_params
#' @inheritParams delay_high
#' @inheritParams dating
#' @inheritDotParams .smooth_LOESS -incidence_input
#' @inheritDotParams .deconvolve_incidence_Richardson_Lucy -incidence_input
#' @inheritDotParams .estimate_Re_EpiEstim_sliding_window -incidence_input
#' @inheritDotParams merge_outputs -output_list -include_index -index_col
#' @inheritDotParams correct_for_partially_observed_data -incidence_data -delay_distribution_final_report
#'
#' @return
#' @export
get_bootstrapped_estimates_from_combined_observations <- function(partially_delayed_incidence,
                                                                  fully_delayed_incidence,
                                                                  smoothing_method = "LOESS",
                                                                  deconvolution_method = "Richardson-Lucy delay distribution",
                                                                  estimation_method = "EpiEstim sliding window",
                                                                  bootstrapping_method = "non-parametric block boostrap",
                                                                  uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                                                  combine_bootstrap_and_estimation_uncertainties = FALSE,
                                                                  N_bootstrap_replicates = 100,
                                                                  delay_until_partial,
                                                                  delay_from_partial_to_full,
                                                                  partial_observation_requires_full_observation = TRUE,
                                                                  ref_date = NULL,
                                                                  time_step = "day",
                                                                  output_Re_only = TRUE,
                                                                  ...) {

  # TODO validate arguments

  # TODO allow for 'partially_delayed_incidence' or 'fully_delayed_incidence' to be NULL,
  # (need to ensure all subsequent functions allow NULL or make if-else)
  # TODO turn get_block_bootstrapped_estimate into a wrapper around this function with partially_delayed_incidence=NULL

  dots_args <- .get_dots_as_list(...)

  index_col <- "idx"
  bootstrap_id_col <- "bootstrap_id"

  # Precompute delay distribution vector or matrix to avoid repeating costly computations needlessly for each bootstrap sample
  delay_distribution_until_partial <- .get_delay_distribution(delay_until_partial,
    n_report_time_steps = .get_input_length(partially_delayed_incidence)
  )
  delay_distribution_partial_to_full <- .get_delay_distribution(delay_from_partial_to_full,
    n_report_time_steps = .get_input_length(fully_delayed_incidence)
  )

  if (partial_observation_requires_full_observation) {
    partially_delayed_incidence <- do.call(
      "correct_for_partially_observed_data",
      c(
        list(
          incidence_data = partially_delayed_incidence,
          delay_distribution_final_report = delay_distribution_partial_to_full
        ),
        .get_shared_args(correct_for_partially_observed_data, dots_args)
      )
    )
  }

  estimate_from_combined_observations_dots_args <- .get_shared_args(
    list(
      .smooth_LOESS,
      .deconvolve_incidence_Richardson_Lucy,
      .estimate_Re_EpiEstim_sliding_window,
      .estimate_Re_EpiEstim_piecewise_constant,
      get_matrix_from_empirical_delay_distr,
      build_delay_distribution
    ),
    dots_args
  )

  original_result <- do.call(
    "estimate_from_combined_observations",
    c(
      list(
        partially_delayed_incidence = partially_delayed_incidence,
        fully_delayed_incidence = fully_delayed_incidence,
        smoothing_method = smoothing_method,
        deconvolution_method = deconvolution_method,
        estimation_method = estimation_method,
        delay_until_partial = delay_distribution_until_partial,
        delay_from_partial_to_full = delay_distribution_partial_to_full,
        partial_observation_requires_full_observation = FALSE,
        ref_date = NULL,
        output_Re_only = FALSE,
        include_index = TRUE,
        index_col = index_col,
        output_HPD = combine_bootstrap_and_estimation_uncertainties
      ),
      estimate_from_combined_observations_dots_args
    )
  )

  if (combine_bootstrap_and_estimation_uncertainties) {
    # Keep the Re estimation HPDs for later
    Re_HPDs <- original_result %>%
      dplyr::select(.data[[index_col]], .data$Re_highHPD, .data$Re_lowHPD)

    # Remove them from the original_result variable
    original_result <- original_result %>%
      dplyr::select(!c(.data$Re_highHPD, .data$Re_lowHPD))
  } else {
    Re_HPDs <- NULL
  }

  original_result[[bootstrap_id_col]] <- 0

  bootstrapping_results <- list(original_result)

  for (i in 1:N_bootstrap_replicates) {
    bootstrapped_partially_delayed_incidence <- do.call(
      "get_bootstrap_replicate",
      c(
        list(
          incidence_data = partially_delayed_incidence,
          bootstrapping_method = bootstrapping_method
        ),
        .get_shared_args(
          list(
            .block_bootstrap,
            .block_bootstrap_overlap_func,
            .smooth_LOESS
          ),
          dots_args
        )
      )
    )

    bootstrapped_fully_delayed_incidence <- do.call(
      "get_bootstrap_replicate",
      c(
        list(
          incidence_data = fully_delayed_incidence,
          bootstrapping_method = bootstrapping_method
        ),
        .get_shared_args(
          list(
            .block_bootstrap,
            .block_bootstrap_overlap_func,
            .smooth_LOESS
          ),
          dots_args
        )
      )
    )

    bootstrapped_estimate <- do.call(
      "estimate_from_combined_observations",
      c(
        list(
          partially_delayed_incidence = bootstrapped_partially_delayed_incidence,
          fully_delayed_incidence = bootstrapped_fully_delayed_incidence,
          smoothing_method = smoothing_method,
          deconvolution_method = deconvolution_method,
          estimation_method = estimation_method,
          delay_until_partial = delay_distribution_until_partial,
          delay_from_partial_to_full = delay_distribution_partial_to_full,
          partial_observation_requires_full_observation = FALSE,
          ref_date = NULL,
          output_Re_only = FALSE,
          include_index = TRUE,
          index_col = index_col,
          outputHPD = FALSE
        ),
        estimate_from_combined_observations_dots_args
      )
    )

    bootstrapped_estimate[[bootstrap_id_col]] <- i

    bootstrapping_results <- c(bootstrapping_results, list(bootstrapped_estimate))
  }

  bootstrapped_estimates <- dplyr::bind_rows(bootstrapping_results)

  original_estimates <- bootstrapped_estimates %>%
    dplyr::filter(.data[[bootstrap_id_col]] == 0)

  bootstrapped_estimates <- bootstrapped_estimates %>%
    dplyr::filter(.data[[bootstrap_id_col]] > 0)

  estimates_with_uncertainty <- do.call(
    "do_uncertainty_summary",
    c(
      list(
        original_values = original_estimates,
        bootstrapped_values = bootstrapped_estimates,
        uncertainty_summary_method = uncertainty_summary_method,
        value_col = "Re_estimate",
        bootstrap_id_col = bootstrap_id_col,
        index_col = index_col,
        output_Re_only = output_Re_only,
        combine_bootstrap_and_estimation_uncertainties = combine_bootstrap_and_estimation_uncertainties,
        Re_HPDs = Re_HPDs
      ),
      .get_shared_args(.summarise_CI_bootstrap, dots_args)
    )
  )

  if (!is.null(ref_date)) {
    estimates_with_uncertainty <- .add_date_column(estimates_with_uncertainty,
      ref_date = ref_date,
      time_step = time_step,
      index_col = index_col,
      keep_index_col = FALSE
    )
  }

  pretty_results <- do.call(
    ".prettify_result",
    c(list(data = estimates_with_uncertainty),
      .get_shared_args(.prettify_result, dots_args))
  )

  return(pretty_results)
}
