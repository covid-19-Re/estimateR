#TODO fill doc
#' Title
#'
#' @param bootstrapped_estimates
#' @param uncertainty_summary_method One of these strings: \itemize{
#' \item{'original estimate - CI from bootstrap estimates'}
#' \item{'bagged mean - CI from bootstrap estimates'}
#' }
#' @param Re_estimate_col string. Name of the column containing Re estimates
#' @param bootstrap_id_col string. Name of the column containing bootstrap samples numbering.
#' Id 0 must correspond to the estimate on the original data.
#' @param time_step string
#'
#' @return dataframe containing Re estimates and confidence interval boundaries.
#' @export
#'
#' @examples
summarise_uncertainty <- function(bootstrapped_estimates,
                                  uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                  Re_estimate_col = "R_mean",
                                  bootstrap_id_col = "bootstrap_id",
                                  time_step = "day"){

  if(uncertainty_summary_method == "original estimate - CI from bootstrap estimates") {
    return(.summarise_CI_bootstrap(bootstrapped_estimates = bootstrapped_estimates,
                                  Re_estimate_col = Re_estimate_col,
                                  bootstrap_id_col = bootstrap_id_col,
                                  time_step = time_step))
  } else if (uncertainty_summary_method == "bagged mean - CI from bootstrap estimates") {
      Re_estimate <- .summarise_bagged_mean(bootstrapped_estimates = bootstrapped_estimates,
                                            include_original_estimates = TRUE,
                                            Re_estimate_col = Re_estimate_col,
                                            bootstrap_id_col = bootstrap_id_col,
                                            time_step = time_step)

      confidence_intervals <- .summarise_CI_bootstrap(bootstrapped_estimates = bootstrapped_estimates,
                                                      Re_estimate_col = Re_estimate_col,
                                                      bootstrap_id_col = bootstrap_id_col,
                                                      time_step = time_step) %>%
                              dplyr::select(-.data$Re_estimate)

      Re_estimate <- Re_estimate %>%
        dplyr::left_join(confidence_intervals, by = "date")

      return(Re_estimate)
  } else {
    stop("Uncertainty summary method is unknown.")
  }
}




#TODO document
#' Title
#'
#' @param bootstrapped_estimates
#' @param alpha
#' @param Re_estimate_col
#' @param bootstrap_id_col
#' @param time_step
#'
#' @return
.summarise_CI_bootstrap <- function(bootstrapped_estimates,
                                   alpha = 0.95,
                                   Re_estimate_col = "R_mean",
                                   bootstrap_id_col = "bootstrap_id",
                                   time_step = "day"){

  high_quantile <- 1-(1-alpha)/2

  original_estimate <- bootstrapped_estimates %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]]), .data[[bootstrap_id_col]] == 0) %>%
    dplyr::select(.data$date, .data[[Re_estimate_col]])

  estimate_with_uncertainty <- bootstrapped_estimates %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]]), !is.na(.data$date), .data[[bootstrap_id_col]] > 0) %>%
    dplyr::select(.data$date, .data[[Re_estimate_col]]) %>%
    dplyr::group_by(.data$date) %>%
    dplyr::summarize(sd_mean = stats::sd(.data[[Re_estimate_col]]),
                     .groups = "drop") %>%
    dplyr::right_join(original_estimate, by = "date") %>%
    dplyr::mutate(CI_down = .data[[Re_estimate_col]] - stats::qnorm(high_quantile)*.data$sd_mean,
                  CI_up = .data[[Re_estimate_col]] + stats::qnorm(high_quantile)*.data$sd_mean) %>%
    dplyr::mutate(CI_down = dplyr::if_else(.data$CI_down < 0, 0, .data$CI_down)) %>%
    dplyr::rename(Re_estimate = .data[[Re_estimate_col]]) %>%
    dplyr::select(-.data$sd_mean) %>%
    tidyr::complete(date = seq.Date(min(.data$date), max(.data$date), by = time_step),
                    fill = list(Re_estimate = NA,
                                CI_down = NA,
                                CI_up = NA))

  return(estimate_with_uncertainty)
}


#TODO document
#' Title
#'
#' @param bootstrapped_estimates
#' @param include_original_estimates
#' @param Re_estimate_col
#' @param bootstrap_id_col
#'
#' @return
.summarise_bagged_mean <- function(bootstrapped_estimates,
                                  include_original_estimates = TRUE,
                                  Re_estimate_col = "R_mean",
                                  bootstrap_id_col = "bootstrap_id",
                                  time_step = "day") {

  if(!include_original_estimates) {
    bootstrapped_estimates <- bootstrapped_estimates %>%
      dplyr::filter(.data[[bootstrap_id_col]] > 0)
  }

  bagged_mean_estimate <- bootstrapped_estimates %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]]), !is.na(.data$date)) %>%
    dplyr::select(.data$date, .data[[Re_estimate_col]]) %>%
    dplyr::group_by(.data$date) %>%
    dplyr::summarize(Re_estimate = mean(.data[[Re_estimate_col]]),
                     .groups = "drop") %>%
    tidyr::complete(date = seq.Date(min(.data$date), max(.data$date), by = time_step),
                    fill = list(Re_estimate = NA))

  return(bagged_mean_estimate)
}

