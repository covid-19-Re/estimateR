#' Convert module output object into tibble
#'
#' @inherit inner_module
#' @param output_name string. Name to be given to the \code{values} column
#' @param index_col string. Name of the index column included in the output.
#' If a \code{ref_date} is passed to the function, the result will contain
#' a \code{date} column instead.
#' Even so, a value must be given to this argument for internal steps.
#' @inherit dating
#'
#' @return tibble
#' @export
make_tibble_from_output <- function(output,
                                    output_name = "value",
                                    index_col = "idx",
                                    ref_date = NULL,
                                    time_step = "day") {
  .are_valid_argument_values(list(
    list(output, "module_input"),
    list(output_name, "string"),
    list(index_col, "string"),
    list(ref_date, "null_or_date"),
    list(time_step, "time_step")
  ))

  tmp_output <- .get_module_input(output)
  indices <- seq(from = .get_offset(tmp_output), by = 1, length.out = length(.get_values(tmp_output)))

  formatted_output <- dplyr::tibble(!!index_col := indices, !!output_name := .get_values(tmp_output))

  if (!is.null(ref_date)) {
    formatted_output <- .add_date_column(
      estimates = formatted_output,
      ref_date = ref_date,
      time_step = time_step,
      index_col = index_col,
      keep_index_col = FALSE
    )
  }

  return(formatted_output)
}


#' Make empty module output object
#'
#' @return empty module output object
.make_empty_module_output <- function() {
  return(list("values" = NA_real_, "index_offset" = 0))
}

#' Transform a result output of a module into a module output 'object'
#'
#' Also takes the module input object \code{input}
#' given to the module to calculate the offset of the output object.
#' The new offset is simply (module input offset) + (shift applied during module operations)
#' The shift applied during module operations is the \code{offset} argument.
#'
#' @param results numeric vector containing output of module operations.
#' @param original_offset integer. Input offset before computations.
#' @param additional_offset integer. Shift resulting from operations performed in module.
#'
#' @inherit module_structure return
.get_module_output <- function(results, original_offset, additional_offset = 0) {
  .are_valid_argument_values(list(
    list(results, "numeric_vector"),
    list(original_offset, "integer"),
    list(additional_offset, "integer")
  ))
  if (length(results) == 0) {
    return(.make_empty_module_output())
  }

  new_offset <- original_offset + additional_offset

  while (is.na(results[1])) {
    results <- results[-1]
    new_offset <- new_offset + 1

    if (length(results) == 0) {
      return(.make_empty_module_output())
    }
  }

  return(list("values" = results, "index_offset" = new_offset))
}

#' Add dates column to dataframe.
#'
#' @param estimates dataframe. Estimates.
#' @param keep_index_col boolean. Keep index column in output?
#' @inherit dating
#' @inherit uncertainty
#'
#' @return estimates dataframe with dates column.
.add_date_column <- function(estimates,
                             ref_date,
                             time_step,
                             index_col = "idx",
                             keep_index_col = FALSE) {
  .are_valid_argument_values(list(
    list(estimates, "estimates", index_col),
    list(ref_date, "date"),
    list(time_step, "time_step"),
    list(index_col, "string"),
    list(keep_index_col, "boolean")
  ))

  dates <- seq.Date(
    from = ref_date + min(estimates[[index_col]]),
    by = time_step,
    along.with = estimates[[index_col]]
  )

  estimates <- estimates %>%
    dplyr::arrange(.data[[index_col]]) %>%
    dplyr::mutate(date = dates)

  estimates <- dplyr::select(estimates, date, tidyselect::everything())

  if (!keep_index_col) {
    estimates <- estimates %>%
      dplyr::select(-.data[[index_col]])
  }

  return(estimates)
}

#' Merge multiple module outputs into tibble
#'
#' Output tibble from list of unsynced module outputs, with an optional date column.
#' The optional \code{ref_date} argument is the starting date of an input with offset 0.
#' In general, this will be the date corresponding to the first entry in the original incidence data.
#' If a reference date is provided with \code{ref_date}, a date column is appended to the tibble,
#' with sequential dates generated with the time step specified by the \code{time_step} parameter.
#'
#'
#' @param output_list named list of module output objects.
#' @param include_index boolean. Include an index column in output?
#' @param index_col string. If \code{include_index} is \code{TRUE},
#' an index column named \code{index_col} is added to the output.
#' @inherit dating
#'
#' @return tibble
#' @export
merge_outputs <- function(output_list,
                          ref_date = NULL,
                          time_step = "day",
                          include_index = is.null(ref_date),
                          index_col = "idx") {
  .are_valid_argument_values(list(
    list(ref_date, "null_or_date"),
    list(time_step, "time_step"),
    list(index_col, "string"),
    list(include_index, "boolean")
  ))
  for (i in 1:length(output_list)) {
    .are_valid_argument_values(list(list(output_list[[i]], "module_input")))
  }

  tibble_list <- lapply(
    1:length(output_list),
    function(i) {
      make_tibble_from_output(
        output = output_list[[i]],
        output_name = names(output_list)[i],
        index_col = index_col
      )
    }
  )

  merged_outputs <- plyr::join_all(tibble_list, by = index_col, type = "full") %>%
    dplyr::arrange(.data[[index_col]])

  if (!is.null(ref_date)) {
    dates <- seq.Date(
      from = ref_date + min(merged_outputs[[index_col]]),
      by = time_step,
      along.with = merged_outputs[[index_col]]
    )
    merged_outputs$date <- dates
    merged_outputs <- dplyr::select(merged_outputs, date, tidyselect::everything())
  }

  if (!include_index) {
    merged_outputs <- dplyr::select(merged_outputs, -.data[[index_col]])
  }

  return(merged_outputs)
}

#' Simplify output object if possible
#'
#' If offset is 0, return only vector containing values.
#' If offset is not zero then \code{output} is returned as is.
#'
#' @param output Module output object
#'
#' @return numeric vector or module output object.
.simplify_output <- function(output) {
  .are_valid_argument_values(list(list(output, "module_input")))

  if (.get_offset(output) == 0) {
    return(.get_values(output))
  } else {
    return(output)
  }
}

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
                             date_col = "date") {
  .are_valid_argument_values(list(
    list(index_col, "string"),
    list(date_col, "string")
  ))

  if (is.data.frame(data)) {
    if (!(index_col %in% names(data) || date_col %in% names(data))) {
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
      dplyr::filter(
        .data[[ref_col]] >= first_row_to_keep,
        .data[[ref_col]] <= last_row_to_keep
      )

    return(cleaned_data)
  } else if (is.list(data)) { # This needs to be checked second, because is.list(dataframe) is TRUE.
    return(.simplify_output(data))
  } else {
    stop("Data must be a list (module output) or dataframe.")
  }
}
