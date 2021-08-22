#' Get the names of the arguments of a function
#'
#' @param func function
#'
#' @return names of the arguments of \code{func}
.get_arg_names <- function(func) {
  return(names(formals(func)))
}

#' Return dots arguments as list.
#'
#' If there are no dots arguments, then return an empty list.
#'
#' @param ... dots arguments.
#'
#' @return list containing dots arguments or empty list
.get_dots_as_list <- function(...) {
  if (...length() > 0) {
    dots_args <- list(...)
  } else {
    dots_args <- list()
  }
  return(dots_args)
}

#' Get the arguments which apply to a function among a given list of arguments.
#'
#' This function is used to find which arguments should be passed
#' to a list of functions among the dot arguments passed to a higher-level function.
#'
#' @param func_list list of functions
#' @param dots_args list of arguments
#'
#' @return elements of \code{dots_args} which can be passed to \code{func_list}
.get_shared_args <- function(func_list, dots_args) {
  if (is.function(func_list)) {
    func_arg_names <- .get_arg_names(func_list)
  } else if (is.list(func_list)) {
    func_arg_names <- unlist(lapply(func_list, .get_arg_names))
  }
  return(dots_args[names(dots_args) %in% func_arg_names])
}
