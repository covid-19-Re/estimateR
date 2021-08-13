#' Transform input data into a module input object
#'
#' The input can be a list containing a \code{values} element
#' and an \code{index_offset} element, potentially among others.
#' It can also be a vector containing numeric values.
#'
#' @param data A module input object or numeric vector.
#'
#' @return module input object.
#' List with a \code{values} and an \code{index_offset} element.
#'
.get_module_input <- function(data) {
  .are_valid_argument_values(list(list(data, "module_input")))

  if (is.list(data)) {
    values <- as.double(data$values)
    index_offset <- data$index_offset
  } else {
    values <- as.double(data)
    index_offset <- 0
  }
  return(list("values" = values, "index_offset" = index_offset))
}

#' Get values from module object
#'
#' This function must be adapted if the module_input_object implementation changes.
#'
#' @param module_object module object.
#'
#' @return vector containing \code{values} only
.get_values <- function(module_object) {
  .are_valid_argument_values(list(list(module_object, "module_input")))
  if (is.list(module_object)) {
    return(module_object$values)
  } else {
    return(module_object)
  }
}


#' Get offset from module object
#'
#' This function must be adapted if the module_input_object implementation changes.
#'
#' @param module_object module object.
#'
#' @return numeric scalar. \code{index_offset} of the \code{module_object}
.get_offset <- function(module_object) {
  .are_valid_argument_values(list(list(module_object, "module_input")))
  if (is.list(module_object)) {
    return(module_object$index_offset)
  } else {
    return(0)
  }
}

#' Get length of values vector in a module input object.
#'
#' @inheritParams inner_module
#'
#' @return integer. length of values vector.
.get_input_length <- function(input) {
  .are_valid_argument_values(list(list(input, "module_input")))
  return(length(.get_values(input)))
}

#' Add values from two module objects
#'
#' Values in the output object are values added from the two objects.
#' The \code{offset} of the output is the maximum offset between the two input objects.
#'
#' @inherit inner_module
#'
#' @export
inner_addition <- function(input_a, input_b) {
  .are_valid_argument_values(list(
    list(input_a, "module_input"),
    list(input_b, "module_input")
  ))

  length_a <- .get_input_length(input_a)
  length_b <- .get_input_length(input_b)

  offset_a <- .get_offset(input_a)
  offset_b <- .get_offset(input_b)

  inner_offset <- max(offset_a, offset_b)
  length_addition <- min(length_a - (inner_offset - offset_a), length_b - (inner_offset - offset_b))

  inner_a <- .get_values(input_a)[seq(from = inner_offset - offset_a + 1, by = 1, length.out = length_addition)]
  inner_b <- .get_values(input_b)[seq(from = inner_offset - offset_b + 1, by = 1, length.out = length_addition)]

  return(.get_module_input(list(values = inner_a + inner_b, index_offset = inner_offset)))
}

#' Add values from two module objects
#'
#' The \code{offset} of the output object is the minimum offset between \code{input_a} and \code{input_b}.
#' @inherit inner_module
#'
#' @export
left_addition <- function(input_a, input_b) {
  .are_valid_argument_values(list(
    list(input_a, "module_input"),
    list(input_b, "module_input")
  ))

  offset_a <- .get_offset(input_a)
  offset_b <- .get_offset(input_b)

  min_offset <- min(offset_a, offset_b)
  padded_input_a <- leftpad_input(input_a, min_offset, padding_value = 0)
  padded_input_b <- leftpad_input(input_b, min_offset, padding_value = 0)

  length_a <- .get_input_length(padded_input_a)
  length_b <- .get_input_length(padded_input_b)

  length_addition <- min(length_a, length_b)

  values_a <- .get_values(padded_input_a)[1:length_addition]
  values_b <- .get_values(padded_input_b)[1:length_addition]

  return(.get_module_input(list(values = values_a + values_b, index_offset = min_offset)))
}

#' Pad values on the left side of input
#'
#' @param new_offset Offset of output.
#' @param padding_value Value to left-pad input with.
#' @inherit inner_module
#'
#' @export
leftpad_input <- function(input, new_offset, padding_value = 0) {
  .are_valid_argument_values(list(
    list(input, "module_input"),
    list(new_offset, "number"),
    list(padding_value, "number")
  ))

  if (new_offset >= .get_offset(input)) {
    return(input)
  } else {
    padded_values <- c(rep(padding_value, length.out = .get_offset(input) - new_offset), .get_values(input))
    return(list(values = padded_values, index_offset = new_offset))
  }
}



