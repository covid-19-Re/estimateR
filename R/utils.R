#TODO make class for input objects (list with values and index_offset fields)

#TODO make a utils that does the input checking (for before transformation and for when it's supposed to have been done)

#TODO fill in doc
#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
get_module_input <- function(data) {
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


#TODO fill in doc
#' Title
#'
#' @param results
#' @param input
#' @param offset
#'
#' @return
get_module_output <- function(results, input, offset = 0) {
  new_offset <- input$index_offset + offset
  if(new_offset == 0){
    # If offset stays 0, output vector with 'results' values
    return(results)
  } else {
    # Else return list with 'values' and 'index_offset
    return(list(values = results, index_offset = new_offset))
  }
}


get_input_length <- function(input) {
  return(length(input$values))
}


#TODO fill in doc
#' Title
#'
#' @param output
#' @param output_name
#'
#' @return
make_tibble_from_output <- function(output, output_name){

  tmp_output <- get_module_input(output)
  indices <- seq(from = tmp_output$index_offset, by = 1, length.out = length(tmp_output$values))

  return(dplyr::tibble(index = indices, !!output_name := tmp_output$values))
}

#Output tibble from list of unsynced module outputs, with an optional date reference
#' Title
#'
#' @param output_list
#' @param ref_date
#' @param time_step
#'
#' @return
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
