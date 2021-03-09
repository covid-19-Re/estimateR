#TODO make class for input objects (list with values and index_offset fields)

#TODO make a utils that does the input checking (for before transformation and for when it's supposed to have been done)

#TODO make a utils function that converts individual module object into a tibble with a date

#TODO clarify the language between module input object and module ouptut object


#' Merge multiple module outputs into tibble
#'
#' Output tibble from list of unsynced module outputs, with an optional date reference.
#' The optional reference date is the starting date of an input with offset 0.
#' In general, this will be the date corresponding to the first entry in the original incidence data.
#' If a reference date is provided, a date column is appended to the tibble,
#' with sequential dates generated with the time step specified by the optional time_step parameter.
#'
#'
#' @param output_list named list of module outputs
#' @param ref_date Date. Reference date.
#' @param time_step string. "day", "2 days", "week", "year"... (see seq.Date help page for details)
#'
#' @return tibble
#' @export
merge_outputs <- function(output_list, ref_date = NULL, time_step = "day"){

  tibble_list <- lapply(1:length(output_list),
                        function(i) {
                          make_tibble_from_output(output = output_list[[i]],
                                                  output_name = names(output_list)[i])
                        })

  merged_outputs<- plyr::join_all(tibble_list, by='index', type='left')

  if( !is.null(ref_date) ) {
    dates <- seq.Date(from = ref_date - min(merged_outputs$index), by = time_step, along.with = merged_outputs$index)
    merged_outputs$date <- dates
  }

  merged_outputs <- dplyr::select(merged_outputs, -index)

  return(merged_outputs)
}


#' Transform input data into a module input object
#'
#' #TODO specify what format data can take
#' @param data
#'
#' @return module input object
.get_module_input <- function(data) {
  #TODO properly check input format
  if (is.list(data)) {
    values <- as.double(data$values)
    index_offset <- data$index_offset
  } else {
    values <- as.double(data)
    index_offset <- 0
  }
  return(list(values = values, index_offset = index_offset))
}


#' Get values from module object
#'
#' This function must be adapted if the module_input_object implementation changes.
#'
#' @param module_input_object
#'
#' @return vector containing values
.get_values <- function(module_object) {
  return(module_object$values)
}


#' Transform a result output of a module into a module output
#'
#' Also takes the module input object to calculate the offset of the output object.
#' The new offset is simply (module input offset) + (shift applied during module operations)
#' The shift applied during module operations is passed as "offset" parameter.
#'
#' @param results numeric vector
#' @param input module input object
#' @param offset integer
#'
#' @return module output object
.get_module_output <- function(results, input, offset = 0) {
  new_offset <- input$index_offset + offset
  if(new_offset == 0){
    # if offset stays 0, output vector with 'results' values
    return(results)
  } else {
    # else return list with 'values' and 'index_offset
    return(list(values = results, index_offset = new_offset))
  }
}


#' Get length of values vector in a module input object.
#'
#' @param input module input object.
#'
#' @return integer. length of values vector.
.get_input_length <- function(input) {
  return(length(input$values))
}


#' Convert module output object into tibble
#'
#' @param output module output object
#' @param output_name string. Name to be given to the values column
#'
#' @return tibble
.make_tibble_from_output <- function(output, output_name){

  tmp_output <- get_module_input(output)
  indices <- seq(from = tmp_output$index_offset, by = 1, length.out = length(tmp_output$values))

  return(dplyr::tibble(index = indices, !!output_name := tmp_output$values))
}



