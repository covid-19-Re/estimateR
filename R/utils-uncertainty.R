#' Summarise the uncertainty obtained from bootstrapping
#'
#' This function takes as input bootstrapped values: Re estimates or others
#' and builds uncertainty intervals from these bootstrapped values.
#'
#' This function is not meant to be needed by regular users.
#' It is exported as it can be useful when building 'pipes':
#' sequences of computations that weave together the different
#' modules provided by \code{estimateR}.
#'
#' @inherit uncertainty
#' @inheritDotParams .summarise_CI_bootstrap
#'
#' @return A dataframe containing Re estimates (column 'Re_estimate')
#' and confidence interval boundaries, with 4 columns like so:
#' \itemize{
#' \item{\code{index_col}, the timestep index column}
#' \item{A column named \code{output_value_col},
#' containing the central values (typically these are Re estimates) }
#' \item{\code{CI_up}, upper limit of the confidence interval}
#' \item{\code{CI_down}, the lower limit of the confidence interval}
#' }
#' 
#' 
#' @export
summarise_uncertainty <- function(bootstrapped_values,
                                  original_values = NULL,
                                  uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
                                  value_col = "Re_estimate",
                                  output_value_col = "Re_estimate",
                                  bootstrap_id_col = "bootstrap_id",
                                  index_col = "idx",
                                  ...) {
  .are_valid_argument_values(list(
    list(bootstrapped_values, "bootstrap_estimates", c(value_col, bootstrap_id_col, index_col)),
    list(uncertainty_summary_method, "uncertainty_summary_method"),
    list(value_col, "string"),
    list(output_value_col, "string"),
    list(bootstrap_id_col, "string"),
    list(index_col, "string")
  ))


  dots_args <- .get_dots_as_list(...)

  bootstrapped_values <- bootstrapped_values %>%
    dplyr::rename(!!output_value_col := .data[[value_col]])

  if (!is.null(original_values)) {
    .are_valid_argument_values(list(list(original_values, "estimates", index_col)))
    original_values <- original_values %>%
      dplyr::rename(!!output_value_col := .data[[value_col]])
  }

  if (uncertainty_summary_method == "original estimate - CI from bootstrap estimates") {
    if (is.null(original_values)) {
      stop("'original_values' must be provided when using uncertainty method
           'original estimate - CI from bootstrap estimates'")
    }

    bootstrap_summary <- do.call(
      ".summarise_CI_bootstrap",
      c(
        list(
          central_values = original_values,
          bootstrapped_values = bootstrapped_values,
          value_col = output_value_col,
          bootstrap_id_col = bootstrap_id_col,
          index_col = index_col
        ),
        .get_shared_args(.summarise_CI_bootstrap, dots_args)
      )
    )
  } else if (uncertainty_summary_method == "bagged mean - CI from bootstrap estimates") {
    central_values <- .summarise_bagged_mean(
      original_values = original_values,
      bootstrapped_values = bootstrapped_values,
      value_col = output_value_col,
      bootstrap_id_col = bootstrap_id_col,
      index_col = index_col
    )

    bootstrap_summary <- do.call(
      ".summarise_CI_bootstrap",
      c(
        list(
          central_values = central_values,
          bootstrapped_values = bootstrapped_values,
          value_col = output_value_col,
          bootstrap_id_col = bootstrap_id_col,
          index_col = index_col
        ),
        .get_shared_args(.summarise_CI_bootstrap, dots_args)
      )
    )
  } else {
    stop("Uncertainty summary method is unknown.")
  }

  return(bootstrap_summary)
}

#' Build a confidence interval from bootstrapped values
#'
#' @inherit uncertainty
#'
#' @return dataframe with 4 columns:
#' \itemize{
#' \item{\code{index_col}, the timestep index column}
#' \item{\code{value_col}, the value input in \code{central_values}}
#' \item{CI_up, the upper limit of the confidence interval}
#' \item{CI_down, the lower limit of the confidence interval}
#' }
.summarise_CI_bootstrap <- function(central_values,
                                    bootstrapped_values,
                                    value_col,
                                    bootstrap_id_col,
                                    index_col,
                                    alpha = 0.95,
                                    prefix_up = "CI_up",
                                    prefix_down = "CI_down") {
  .are_valid_argument_values(list(
    list(central_values, "estimates", index_col),
    list(bootstrapped_values, "bootstrap_estimates", c(value_col, bootstrap_id_col, index_col)),
    list(value_col, "string"),
    list(bootstrap_id_col, "string"),
    list(index_col, "string"),
    list(alpha, "numeric_between_zero_one"),
    list(prefix_up, "string"),
    list(prefix_down, "string")
  ))

  CI_down <- paste(prefix_down, value_col, sep = "_")
  CI_up <- paste(prefix_up, value_col, sep = "_")

  high_quantile <- 1 - (1 - alpha) / 2

  central_values <- central_values %>%
    dplyr::select(.data[[index_col]], .data[[value_col]]) %>%
    dplyr::filter(!is.na(.data[[value_col]]))

  value_with_uncertainty <- bootstrapped_values %>%
    dplyr::select(.data[[index_col]], .data[[value_col]]) %>%
    dplyr::filter(!is.na(.data[[value_col]])) %>%
    dplyr::group_by(.data[[index_col]]) %>%
    dplyr::summarize(
      sd_mean = stats::sd(.data[[value_col]]),
      .groups = "drop"
    ) %>%
    dplyr::right_join(central_values, by = index_col) %>%
    dplyr::mutate(
      !!CI_down := .data[[value_col]] - stats::qnorm(high_quantile) * .data$sd_mean,
      !!CI_up := .data[[value_col]] + stats::qnorm(high_quantile) * .data$sd_mean
    ) %>%
    dplyr::mutate(!!CI_down := dplyr::if_else(.data[[CI_down]] < 0, 0, .data[[CI_down]])) %>%
    dplyr::select(-.data$sd_mean) %>%
    tidyr::complete(!!index_col := seq(min(.data[[index_col]]), max(.data[[index_col]])))

  return(value_with_uncertainty)
}

#' Compute bagged mean from bootstrapped replicates
#'
#' If \code{original_values} are included,
#' these values are included in the mean computation
#' along with the \code{bootstrapped_values}.
#'
#' @inherit uncertainty
#'
#' @return a dataframe containing a time step index column named \code{index_col}
#' and a column containing bagged mean values called \code{value_col}
.summarise_bagged_mean <- function(bootstrapped_values,
                                   original_values = NULL,
                                   value_col,
                                   bootstrap_id_col,
                                   index_col) {
  .are_valid_argument_values(list(
    list(bootstrapped_values, "bootstrap_estimates", c(value_col, bootstrap_id_col, index_col)),
    list(value_col, "string"),
    list(bootstrap_id_col, "string"),
    list(index_col, "string")
  ))


  bootstrapped_values <- bootstrapped_values %>%
    dplyr::select(.data[[index_col]], .data[[value_col]])

  if (!is.null(original_values)) {
    .are_valid_argument_values(list(list(original_values, "estimates", index_col)))

    original_values <- original_values %>%
      dplyr::select(.data[[index_col]], .data[[value_col]])

    bootstrapped_values <- bootstrapped_values %>%
      dplyr::bind_rows(original_values)
  }

  bagged_mean_value <- bootstrapped_values %>%
    dplyr::filter(!is.na(.data[[value_col]])) %>%
    dplyr::group_by(.data[[index_col]]) %>%
    dplyr::summarize(!!value_col := mean(.data[[value_col]]),
      .groups = "drop"
    ) %>%
    tidyr::complete(!!index_col := seq(min(.data[[index_col]]), max(.data[[index_col]])))

  return(bagged_mean_value)
}


#' Utility for summarising uncertainty
#'
#' This function is not meant to be used by typical users.
#' It can be used to build custom pipes with \code{estimateR}.
#'
#' @inherit uncertainty
#' @inherit pipe_params
#' @inheritDotParams .summarise_CI_bootstrap
#'
#' @return A dataframe containing Re estimates (column 'Re_estimate')
#' and confidence interval boundaries, with 4 columns like so:
#' \itemize{
#' \item{\code{index_col}, the timestep index column}
#' \item{A column named \code{output_value_col},
#' containing the central values (typically these are Re estimates) }
#' \item{\code{CI_up}, upper limit of the confidence interval}
#' \item{\code{CI_down}, the lower limit of the confidence interval}
#' }
do_uncertainty_summary <- function(original_values,
                                   bootstrapped_values,
                                   uncertainty_summary_method,
                                   value_col,
                                   bootstrap_id_col,
                                   index_col,
                                   output_Re_only,
                                   combine_bootstrap_and_estimation_uncertainties = FALSE,
                                   Re_HPDs = NULL,
                                   ...) {
  dots_args <- .get_dots_as_list(...)

  CI_down_col_name <- paste0("CI_down_", value_col)
  CI_up_col_name <- paste0("CI_up_", value_col)

  if (output_Re_only) {
    estimates_with_uncertainty <- do.call(
      "summarise_uncertainty",
      c(
        list(
          original_values = original_values,
          bootstrapped_values = bootstrapped_values,
          uncertainty_summary_method = uncertainty_summary_method,
          value_col = value_col,
          output_value_col = value_col,
          bootstrap_id_col = bootstrap_id_col,
          index_col = index_col
        ),
        .get_shared_args(.summarise_CI_bootstrap, dots_args)
      )
    )

    CI_down_col_name <- paste0("CI_down_", value_col)
    CI_up_col_name <- paste0("CI_up_", value_col)

    if (combine_bootstrap_and_estimation_uncertainties) {
      estimates_with_uncertainty <- dplyr::full_join(estimates_with_uncertainty,
        Re_HPDs,
        by = index_col
      ) %>%
        dplyr::mutate(
          !!CI_down_col_name := dplyr::if_else(.data[[CI_down_col_name]] > .data$Re_lowHPD,
            .data$Re_lowHPD, .data[[CI_down_col_name]]
          ),
          !!CI_up_col_name := dplyr::if_else(.data[[CI_up_col_name]] < .data$Re_highHPD,
            .data$Re_highHPD, .data[[CI_up_col_name]]
          )
        ) %>%
        dplyr::select(!c(.data$Re_lowHPD, .data$Re_highHPD))
    }
  } else {
    cols_to_summarise <- names(bootstrapped_values)
    cols_to_summarise <- cols_to_summarise[!cols_to_summarise %in% c(index_col, bootstrap_id_col)]

    summaries <- lapply(cols_to_summarise, function(col_x) {
      bootstrapped_estimates_of_interest <- bootstrapped_values %>%
        dplyr::select(.data[[col_x]], .data[[index_col]], .data[[bootstrap_id_col]])

      original_estimates_of_interest <- original_values %>%
        dplyr::select(.data[[col_x]], .data[[index_col]], .data[[bootstrap_id_col]])

      do.call(
        "summarise_uncertainty",
        c(
          list(
            original_values = original_estimates_of_interest,
            bootstrapped_values = bootstrapped_estimates_of_interest,
            uncertainty_summary_method = uncertainty_summary_method,
            value_col = col_x,
            output_value_col = col_x,
            bootstrap_id_col = bootstrap_id_col,
            index_col = index_col
          ),
          .get_shared_args(.summarise_CI_bootstrap, dots_args)
        )
      )
    })

    estimates_with_uncertainty <- summaries %>%
      purrr::reduce(dplyr::full_join, by = index_col)

    if (combine_bootstrap_and_estimation_uncertainties) {
      bootstrapped_CI_down_col_name <- paste0("bootstrapped_CI_down_", value_col)
      bootstrapped_CI_up_col_name <- paste0("bootstrapped_CI_up_", value_col)


      estimates_with_uncertainty <- dplyr::full_join(estimates_with_uncertainty, Re_HPDs,
        by = index_col
      ) %>%
        dplyr::mutate(
          !!bootstrapped_CI_down_col_name := .data[[CI_down_col_name]],
          !!bootstrapped_CI_up_col_name := .data[[CI_up_col_name]]
        ) %>%
        dplyr::mutate(
          !!CI_down_col_name := dplyr::if_else(.data[[CI_down_col_name]] > .data$Re_lowHPD,
            .data$Re_lowHPD, .data[[CI_down_col_name]]
          ),
          !!CI_up_col_name := dplyr::if_else(.data[[CI_up_col_name]] < .data$Re_highHPD,
            .data$Re_highHPD, .data[[CI_up_col_name]]
          )
        )
    }
  }

  return(estimates_with_uncertainty)
}
