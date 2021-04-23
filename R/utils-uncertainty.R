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
                                  original_estimates = NA,
                                  uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                  Re_estimate_col = "R_mean",
                                  bootstrap_id_col = "bootstrap_id",
                                  time_step = "day"){

  if(Re_estimate_col %!in% names(bootstrapped_estimates)) {
    stop(paste0("Missing ", Re_estimate_col, " column in 'bootstrapped estimates' argument,
                or 'Re_estimate_col' was not set to the corresponding column name."))
  }

  if(bootstrap_id_col %!in% names(bootstrapped_estimates)) {
    stop(paste0("Missing ", bootstrap_id_col, " column in 'bootstrapped estimates' argument,
                or 'bootstrap_id_col' was not set to the corresponding column name."))
  }

  old_Re_estimate_col <- Re_estimate_col
  Re_estimate_col <- "Re_estimate"

  bootstrapped_estimates <- bootstrapped_estimates %>%
    dplyr::rename(!!Re_estimate_col := .data[[old_Re_estimate_col]],
           bootstrap_id = .data[[bootstrap_id_col]])

  if(!is.na(original_estimates)) {
    original_estimates <- original_estimates %>%
      dplyr::rename(!!Re_estimate_col := .data[[old_Re_estimate_col]])
  }

  if(uncertainty_summary_method == "original estimate - CI from bootstrap estimates") {

    if(is.na(original_estimates)) {
      stop("'original_estimates' must be provided when using uncertainty method
           'original estimate - CI from bootstrap estimates'")
    }

    original_estimates <- original_estimates %>%
      dplyr::rename(!!Re_estimate_col := .data[[Re_estimate_col]])

    return(.summarise_CI_bootstrap(central_estimates = original_estimates,
                                  bootstrapped_estimates = bootstrapped_estimates,
                                  Re_estimate_col = Re_estimate_col,
                                  bootstrap_id_col = bootstrap_id_col,
                                  time_step = time_step))

  } else if (uncertainty_summary_method == "bagged mean - CI from bootstrap estimates") {

      central_estimates <- .summarise_bagged_mean(original_estimates = original_estimates,
                                            bootstrapped_estimates = bootstrapped_estimates,
                                            Re_estimate_col = Re_estimate_col,
                                            bootstrap_id_col = bootstrap_id_col,
                                            time_step = time_step)

      Re_estimate <- .summarise_CI_bootstrap(central_estimates = central_estimates,
                                             bootstrapped_estimates = bootstrapped_estimates,
                                             Re_estimate_col = Re_estimate_col,
                                             bootstrap_id_col = bootstrap_id_col,
                                             time_step = time_step)

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
.summarise_CI_bootstrap <- function(central_estimates,
                                    bootstrapped_estimates,
                                   alpha = 0.95,
                                   Re_estimate_col = Re_estimate_col,
                                   bootstrap_id_col = bootstrap_id_col,
                                   time_step = "day"){

  high_quantile <- 1-(1-alpha)/2

  central_estimates <- central_estimates %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]]))

  estimate_with_uncertainty <- bootstrapped_estimates %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]]), !is.na(.data$date)) %>%
    dplyr::select(.data$date, .data[[Re_estimate_col]]) %>%
    dplyr::group_by(.data$date) %>%
    dplyr::summarize(sd_mean = stats::sd(.data[[Re_estimate_col]]),
                     .groups = "drop") %>%
    dplyr::right_join(central_estimates, by = "date") %>%
    dplyr::mutate(CI_down = .data[[Re_estimate_col]] - stats::qnorm(high_quantile)*.data$sd_mean,
                  CI_up = .data[[Re_estimate_col]] + stats::qnorm(high_quantile)*.data$sd_mean) %>%
    dplyr::mutate(CI_down = dplyr::if_else(.data$CI_down < 0, 0, .data$CI_down)) %>%
    dplyr::select(-.data$sd_mean) %>%
    tidyr::complete(date = seq.Date(min(.data$date), max(.data$date), by = time_step))

  return(estimate_with_uncertainty)
}


#TODO document
#' Title
#'
#' @param bootstrapped_estimates
#' @param Re_estimate_col
#' @param bootstrap_id_col
#'
#' @return
.summarise_bagged_mean <- function(original_estimates,
                                  bootstrapped_estimates,
                                  Re_estimate_col = Re_estimate_col,
                                  bootstrap_id_col = bootstrap_id_col,
                                  time_step = "day") {


  bootstrapped_estimates <- bootstrapped_estimates %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]]), !is.na(.data$date)) %>%
    dplyr::select(.data$date, .data[[Re_estimate_col]])

  if(!is.na(original_estimates)) {
      original_estimates <- original_estimates %>%
        dplyr::filter(!is.na(.data[[Re_estimate_col]]), !is.na(.data$date)) %>%
        dplyr::select(.data$date, .data[[Re_estimate_col]])

      bootstrapped_estimates <- bootstrapped_estimates %>%
        dplyr::bind_rows(original_estimates)
  }

  bagged_mean_estimate <- bootstrapped_estimates %>%
    dplyr::group_by(.data$date) %>%
    dplyr::summarize(!!Re_estimate_col := mean(.data[[Re_estimate_col]]),
                     .groups = "drop") %>%
    tidyr::complete(date = seq.Date(min(.data$date), max(.data$date), by = time_step),
                    fill = list(Re_estimate = NA))

  return(bagged_mean_estimate)
}

