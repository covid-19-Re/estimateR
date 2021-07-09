#TODO doc
do_uncertainty_summary <- function(original_values,
                                   bootstrapped_values,
                                   uncertainty_summary_method,
                                   value_col,
                                   bootstrap_id_col,
                                   index_col,
                                   output_Re_only,
                                   combine_bootstrap_and_estimation_uncertainties = FALSE,
                                   Re_HPDs = NULL,
                                   ...){

  dots_args <- .get_dots_as_list(...)

  CI_down_col_name <- paste0("CI_down_", value_col)
  CI_up_col_name <- paste0("CI_up_", value_col)

  if(output_Re_only) {
#TODO remove leading NAs from output
    estimates_with_uncertainty <- do.call(
      'summarise_uncertainty',
      c(list(original_values = original_values,
             bootstrapped_values = bootstrapped_values,
             uncertainty_summary_method = uncertainty_summary_method,
             value_col = value_col,
             output_value_col = value_col,
             bootstrap_id_col = bootstrap_id_col,
             index_col = index_col),
        .get_shared_args(.summarise_CI_bootstrap, dots_args))
      )

    CI_down_col_name <- paste0("CI_down_", value_col)
    CI_up_col_name <- paste0("CI_up_", value_col)

    if(combine_bootstrap_and_estimation_uncertainties) {
      estimates_with_uncertainty <- dplyr::full_join(estimates_with_uncertainty,
                                                     Re_HPDs,
                                                     by = index_col) %>%
        dplyr::mutate(!!CI_down_col_name := dplyr::if_else(.data[[CI_down_col_name]] > .data$Re_lowHPD,
                                                           .data$Re_lowHPD, .data[[CI_down_col_name]]),
                      !!CI_up_col_name := dplyr::if_else(.data[[CI_up_col_name]] < .data$Re_highHPD,
                                                         .data$Re_highHPD, .data[[CI_up_col_name]])) %>%
        dplyr::select(!c(.data$Re_lowHPD, .data$Re_highHPD))
    }

  } else {

    cols_to_summarise <- names(bootstrapped_values)
    cols_to_summarise <- cols_to_summarise[!cols_to_summarise %in% c(index_col, bootstrap_id_col) ]

    summaries <- lapply(cols_to_summarise, function(col_x){
      bootstrapped_estimates_of_interest <- bootstrapped_values %>%
        dplyr::select(.data[[col_x]], .data[[index_col]], .data[[bootstrap_id_col]])

      original_estimates_of_interest <- original_values %>%
        dplyr::select(.data[[col_x]], .data[[index_col]], .data[[bootstrap_id_col]])

      do.call( 'summarise_uncertainty',
        c(list(original_values = original_estimates_of_interest,
               bootstrapped_values = bootstrapped_estimates_of_interest,
               uncertainty_summary_method = uncertainty_summary_method,
               value_col = col_x,
               output_value_col = col_x,
               bootstrap_id_col = bootstrap_id_col,
               index_col = index_col),
          .get_shared_args(.summarise_CI_bootstrap, dots_args))
      )
    })

    estimates_with_uncertainty <- summaries %>%
      purrr::reduce(dplyr::full_join, by = index_col)

    if(combine_bootstrap_and_estimation_uncertainties) {

      bootstrapped_CI_down_col_name <- paste0("bootstrapped_CI_down_", value_col)
      bootstrapped_CI_up_col_name <- paste0("bootstrapped_CI_up_", value_col)


      estimates_with_uncertainty <- dplyr::full_join(estimates_with_uncertainty, Re_HPDs,
                                                     by = index_col) %>%
        dplyr::mutate(!!bootstrapped_CI_down_col_name := .data[[CI_down_col_name]],
                      !!bootstrapped_CI_up_col_name := .data[[CI_up_col_name]]) %>%
        dplyr::mutate(CI_down_Re_estimate = dplyr::if_else(.data$CI_down_Re_estimate > .data$Re_lowHPD,
                                                           .data$Re_lowHPD, .data$CI_down_Re_estimate),
                      CI_up_Re_estimate = dplyr::if_else(.data$CI_up_Re_estimate < .data$Re_highHPD,
                                                         .data$Re_highHPD, .data$CI_up_Re_estimate))
    }
  }

  return(estimates_with_uncertainty)
}

