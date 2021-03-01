
#TODO fill in doc
#' Title
#'
#' @param incidence_data
#' @param estimation_method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
estimate_Re <- function( incidence_data, estimation_method = "EpiEstim sliding window", ... ) {

  input <- get_module_input(incidence_data)

  Re_estimate <- dplyr::case_when(
    estimation_method == "EpiEstim sliding window" ~ estimate_Re_EpiEstim_sliding_window(input, ... ),
    TRUE ~ rep(NA_real_, length.out = get_input_length(input))
  )
  return(Re_estimate)
}


#TODO fill in doc
#' Title
#'
#' @param incidence_input
#' @param minimul_cumul_incidence Numeric scalar. Minimum number of cumulated infections before starting the Re estimation
#' @param estimation_window
#' @param mean_serial_interval
#' @param std_serial_interval
#' @param mean_Re_prior
#'
#' @return
estimate_Re_EpiEstim_sliding_window <- function(incidence_input,
                        minimul_cumul_incidence = 0,
                        estimation_window = 3,
                        mean_serial_interval = 4.8,
                        std_serial_interval  = 2.3,
                        mean_Re_prior = 1) {

  incidence_vector <- incidence_input$values

  # Ensure that incidence_vector has no NA or negative values
  #TODO make these checks in the get_module_input function
  incidence_vector <- dplyr::if_else(is.na(incidence_vector) | incidence_vector < 0, 0, incidence_vector)

  offset <- which(cumsum(incidence_vector) >= minimul_cumul_incidence)[1]
  # offset needs to be at least two for EpiEstim
  offset <- max(2, offset)

  right_bound <- length(incidence_vector) - (estimation_window - 1)

  if (is.na(offset) | right_bound < offset) {
    # No valid data point, return empty estimate
    return(c())
  }

  # Computation intervals corresponding to every position of the sliding window
  t_start <- seq(offset, right_bound)
  t_end <- t_start + estimation_window - 1

  R_instantaneous <- EpiEstim::estimate_R(
    incidence_vector,
    method = "parametric_si",
    config = EpiEstim::make_config(
      list(
        mean_si = mean_serial_interval,
        std_si = std_serial_interval,
        t_start = t_start,
        t_end = t_end,
        mean_prior = mean_Re_prior)
    )
  )


  additional_offset <- t_end[1] - 1
  R_mean <- get_module_output(R_instantaneous$R$`Mean(R)`,
                              input,
                              additional_offset)

  return(R_mean)
}
