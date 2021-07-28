#' Prettify results of pipe functions by removing leading and tailing NAs
#'
#' @param data Module object or dataframe.
#' If dataframe, must contain a column named \code{index_col}
#' or \code{date_col} (or both),
#' @param index_col string. Name of the index column.
#' @param date_col string. Name of the date column.
#'
#' @return The input dataframe without leading NA rows.
.prettify_result <- function(data,
                             index_col = "idx",
                             date_col = "date"){
  .are_valid_argument_values(list(
    list(index_col, "string"),
    list(date_col, "string")
  ))

  if(is.data.frame(data)){
    if(!(index_col %in% names(data) || date_col %in% names(data))){
      stop("data argument must contain an index column or date column (or both).")
    }

    ref_col <- ifelse(date_col %in% names(data), date_col, index_col)

    # Remove leading rows with only NAs
    first_row_to_keep <- data %>%
      dplyr::arrange(.data[[ref_col]]) %>%
      dplyr::filter(dplyr::if_any(!dplyr::any_of(c(index_col, date_col)), ~ !is.na(.))) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::pull(.data[[ref_col]])

    last_row_to_keep <- data %>%
      dplyr::arrange(dplyr::desc(.data[[ref_col]])) %>%
      dplyr::filter(dplyr::if_any(!dplyr::any_of(c(index_col, date_col)), ~ !is.na(.))) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::pull(.data[[ref_col]])


    cleaned_data <- data %>%
      dplyr::filter(.data[[ref_col]] >= first_row_to_keep,
                    .data[[ref_col]] <= last_row_to_keep)

    return(cleaned_data)
  } else if(is.list(data)) { # This needs to be checked second, because is.list(dataframe) is TRUE.
    return(.simplify_output(data))
  } else {
    stop("Data must be a list (module output) or dataframe.")
  }
}
