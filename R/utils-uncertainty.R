#' Summarise the uncertainty obtained from bootstrapping
#'
#' @param uncertainty_summary_method One of these options:
#' \itemize{
#' \item{'original estimate - CI from bootstrap estimates'.
#' The confidence interval is built using bootstrapped estimates
#' and centered around the original estimates.}
#' \item{'bagged mean - CI from bootstrap estimates'.
#' The confidence interval is built using bootstrapped estimates
#' and centered around the mean of bootstrapped estimates and original estimates.}
#' }
#' @inherit uncertainty
#' @inheritDotParams .summarise_CI_bootstrap
#'
#' @return dataframe containing Re estimates (column 'Re_estimate')
#' and confidence interval boundaries, with 4 columns like so:
#' \itemize{
#' \item{\code{index_col}, the timestep index column}
#' \item{\code{Re_estimate_col}, Re estimates}
#' \item{CI_up, the upper limit of the confidence interval}
#' \item{CI_down, the lower limit of the confidence interval}
#' }
#' @export
summarise_uncertainty <- function(bootstrapped_estimates,
                                  original_estimates = NULL,
                                  uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                  Re_estimate_col = "R_mean",
                                  bootstrap_id_col = "bootstrap_id",
                                  index_col = "idx",
                                  ...){

  if(Re_estimate_col %!in% names(bootstrapped_estimates)) {
    stop(paste0("Missing ", Re_estimate_col, " column in 'bootstrapped estimates' argument,
                or 'Re_estimate_col' was not set to the corresponding column name."))
  }

  if(bootstrap_id_col %!in% names(bootstrapped_estimates)) {
    stop(paste0("Missing ", bootstrap_id_col, " column in 'bootstrapped estimates' argument,
                or 'bootstrap_id_col' was not set to the corresponding column name."))
  }

  if(index_col %!in% names(bootstrapped_estimates)) {
    stop(paste0("Missing ", index_col, " column in 'bootstrapped estimates' argument,
                or 'index_col' was not set to the corresponding column name."))
  }

  dots_args <- .get_dots_as_list(...)

  old_Re_estimate_col <- Re_estimate_col
  Re_estimate_col <- "Re_estimate"

  bootstrapped_estimates <- bootstrapped_estimates %>%
    dplyr::rename(!!Re_estimate_col := .data[[old_Re_estimate_col]])

  if(!is.null(original_estimates)) {
    original_estimates <- original_estimates %>%
      dplyr::rename(!!Re_estimate_col := .data[[old_Re_estimate_col]])
  }

  if(uncertainty_summary_method == "original estimate - CI from bootstrap estimates") {

    if( is.null(original_estimates) ) {
      stop("'original_estimates' must be provided when using uncertainty method
           'original estimate - CI from bootstrap estimates'")
    }

    CI_bootstrap <- do.call(
      '.summarise_CI_bootstrap',
      c(list(central_estimates = original_estimates,
             bootstrapped_estimates = bootstrapped_estimates,
             Re_estimate_col = Re_estimate_col,
             bootstrap_id_col = bootstrap_id_col,
             index_col = index_col),
        .get_shared_args(.summarise_CI_bootstrap, dots_args))
    )

    return(CI_bootstrap)

  } else if (uncertainty_summary_method == "bagged mean - CI from bootstrap estimates") {

    central_estimates <- .summarise_bagged_mean(original_estimates = original_estimates,
                                                bootstrapped_estimates = bootstrapped_estimates,
                                                Re_estimate_col = Re_estimate_col,
                                                bootstrap_id_col = bootstrap_id_col,
                                                index_col = index_col)

    Re_estimate <- do.call(
      '.summarise_CI_bootstrap',
      c(list(central_estimates = central_estimates,
             bootstrapped_estimates = bootstrapped_estimates,
             Re_estimate_col = Re_estimate_col,
             bootstrap_id_col = bootstrap_id_col,
             index_col = index_col),
        .get_shared_args(.summarise_CI_bootstrap, dots_args))
    )

    return(Re_estimate)
  } else {
    stop("Uncertainty summary method is unknown.")
  }
}

#' Build a confidence interval from bootstrapped estimates
#'
#' @inherit uncertainty
#'
#' @return dataframe with 4 columns:
#' \itemize{
#' \item{\code{index_col}, the timestep index column}
#' \item{\code{Re_estimate_col}, the estimate input in \code{central_estimates}}
#' \item{CI_up, the upper limit of the confidence interval}
#' \item{CI_down, the lower limit of the confidence interval}
#' }
.summarise_CI_bootstrap <- function(central_estimates,
                                    bootstrapped_estimates,
                                    Re_estimate_col = Re_estimate_col,
                                    bootstrap_id_col = bootstrap_id_col,
                                    index_col = index_col,
                                    alpha = 0.95){

  #TODO proper validation of input (check that numeric between 0 and 1,
  # strings and dataframes with the right columns and with no NA in index_col)

  if(any(is.na(bootstrapped_estimates[[index_col]]))) {
    stop(paste("NA value(s) in column", index_col, "in", deparse(substitute(bootstrapped_estimates))))
  }

  if(any(is.na(central_estimates[[index_col]]))) {
    stop(paste("NA value(s) in column", index_col, "in", deparse(substitute(central_estimates))))
  }

  high_quantile <- 1-(1-alpha)/2

  central_estimates <- central_estimates %>%
    dplyr::select(.data[[index_col]], .data[[Re_estimate_col]]) %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]]))

  estimate_with_uncertainty <- bootstrapped_estimates %>%
    dplyr::select(.data[[index_col]], .data[[Re_estimate_col]]) %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]])) %>%
    dplyr::group_by(.data[[index_col]]) %>%
    dplyr::summarize(sd_mean = stats::sd(.data[[Re_estimate_col]]),
                     .groups = "drop") %>%
    dplyr::right_join(central_estimates, by = index_col) %>%
    dplyr::mutate(CI_down = .data[[Re_estimate_col]] - stats::qnorm(high_quantile)*.data$sd_mean,
                  CI_up = .data[[Re_estimate_col]] + stats::qnorm(high_quantile)*.data$sd_mean) %>%
    dplyr::mutate(CI_down = dplyr::if_else(.data$CI_down < 0, 0, .data$CI_down)) %>%
    dplyr::select(-.data$sd_mean) %>%
    tidyr::complete(!!index_col := seq(min(.data[[index_col]]), max(.data[[index_col]])))

  return(estimate_with_uncertainty)
}

#' Compute bagged mean from bootstrapped replicates
#'
#' If \code{original_estimates} are included,
#' these estimates are included in the mean computation
#' along with the \code{bootstrapped_estimates}.
#'
#' @inherit uncertainty
#'
#' @return a dataframe containing an timestep index column named \code{index_col}
#' and a column containing bagged mean estimates called \code{Re_estimate_col}
.summarise_bagged_mean <- function(bootstrapped_estimates,
                                   original_estimates = NULL,
                                   Re_estimate_col = Re_estimate_col,
                                   bootstrap_id_col = bootstrap_id_col,
                                   index_col = index_col) {

  #TODO proper validation of input (check that strings and dataframes with the right columns and with no NA in index_col)

  if(any(is.na(bootstrapped_estimates[[index_col]]))) {
    stop(paste("NA value(s) in column", index_col, "in", deparse(substitute(bootstrapped_estimates))))
  }

  bootstrapped_estimates <- bootstrapped_estimates %>%
    dplyr::select(.data[[index_col]], .data[[Re_estimate_col]])

  if(!is.null(original_estimates)) {
    if(any(is.na(original_estimates[[index_col]]))) {
      stop(paste("NA value(s) in column", index_col, "in", deparse(substitute(original_estimates))))
    }

    original_estimates <- original_estimates %>%
      dplyr::select(.data[[index_col]], .data[[Re_estimate_col]])

    bootstrapped_estimates <- bootstrapped_estimates %>%
      dplyr::bind_rows(original_estimates)
  }

  bagged_mean_estimate <- bootstrapped_estimates %>%
    dplyr::filter(!is.na(.data[[Re_estimate_col]])) %>%
    dplyr::group_by(.data[[index_col]]) %>%
    dplyr::summarize(!!Re_estimate_col := mean(.data[[Re_estimate_col]]),
                     .groups = "drop") %>%
    tidyr::complete(!!index_col := seq(min(.data[[index_col]]), max(.data[[index_col]])))

  return(bagged_mean_estimate)
}
