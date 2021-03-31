#TODO write deconvolution step with the convolution of a matrix with vector. That should apply to Swiss non-onset data
#TODO write deconvolution step that takes into account that data is onset data (when it is). That should apply to Swiss onset data

#TODO add a function to deal with empirical delay distribution and build delay distribution matrix
## make sure it can deal with Spanish data specificity

#TODO transform make_ecdf_from_two_gammas into a more general function summing draws from n distributions (allow type of distribution to be changed)

#TODO make make_ecdf_from_empirical_data_and_gamma more general, draws can be from different types of distribution.
## Check why gamma_draws is an input and not drawn inside the function
## Check why not Vectorized ecdf output as in make_ecdf_from_two_gammas

#TODO think about whether utilities like '.get_matrix_from_single_delay_distr' need to be exported


#' Make square delay distribution matrix from vector of delay distributions.
#'
#' @param delay_distribution numeric vector
#' @param N integer. Dimension of output matrix
#'
#' @return square numeric matrix
.get_matrix_from_single_delay_distr <- function(delay_distribution, N) {

  if(N >= length(delay_distribution)) {
    delay_distribution <- c(delay_distribution, rep(0, times = N - length(delay_distribution)))
  }

  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), delay_distribution[1:(N - i + 1)])
  }

  return(delay_distribution_matrix)
}

#TODO fill documentation
#TODO test
#maybe merge with .get_matrix_from_single_delay_distr by adding N parm and checking if list or unique vector
#' Build delay distribution matrix from list of delay distribution vectors
#'
#' @param waiting_time_distribution_list
#'
#' @return
.get_matrix_from_waiting_time_distributions <- function(waiting_time_distribution_list){
  N <- length(waiting_time_distribution_list)
  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)

  for(i in 1:N){

    waiting_time_distr <- waiting_time_distribution_list[[i]]
    if(length(waiting_time_distr) < N - i + 1) {
      waiting_time_distr <- c(waiting_time_distr, rep(0, times = N - i + 1 - length(waiting_time_distr)))
    }
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), waiting_time_distr[1:(N - i + 1)])
  }

  return(delay_distribution_matrix)
}

#TODO add details on the discretization
#TODO fill documentation
#TODO test (test that vector sums up to 1)
#' Build a delay distribution vector
#'
#' Only allows for gamma distributions for now.
#' @param parm1 numeric. If \code{distribution_type=="gamma"}, \code{parm1} is the shape parameter.
#' @param parm2 numeric. If \code{distribution_type=="gamma"}, \code{parm2} is the scale parameter.
#' @param distribution_type string. Options are "gamma".
#' @param max_quantile numeric value between 0 and 1. TODO write what max_quantile does
#'
#' @return numeric vector.
#' @export
#'
#' @examples
#' #TODO add example
build_delay_distribution <- function(parm1,
                                     parm2,
                                     distribution_type = "gamma",
                                     max_quantile = 0.999){
  if(distribution_type == "gamma") {
    # Take the right boundary of the delay distribution vector
    right_boundary <- ceiling(stats::qgamma(max_quantile, shape = parm1, scale = parm2)) + 1
    right_boundary <- max(right_boundary, 2) # Set the right boundary to at least two

    cdf_values <- stats::pgamma(c(0, seq(from = 0.5, to = right_boundary, by = 1)), shape = parm1, scale = parm2)
    distribution_vector <- diff(cdf_values)
  } else {
    #TODO throw error
    return(NA)
  }

  return(distribution_vector)
}


#TODO fill documentation
#' Title
#'
#' @param vector_a numeric
#' @param vector_b numeric
#'
#' @return numeric vector
.convolve_delay_distribution_vectors <- function(vector_a, vector_b){
  # Right-pad vectors with zeroes to bring them to the same length

  final_length <- length(vector_b) + length(vector_a)

  vector_a <- c(vector_a, rep(0,times = final_length - length(vector_a)))
  vector_b <- c(vector_b, rep(0,times = final_length - length(vector_b)))

  vector_c <- rep(0, times = final_length)

  for(i in 1 : final_length) {
    reversed_vector_b <- rev(vector_b[1:i]) # Reverse vector_b truncated at index i
    reversed_vector_b <- c(reversed_vector_b, rep(0, times = final_length - i)) # Right-pad with zeroes
    vector_c[i] <- vector_a %*% reversed_vector_b # Compute dot product between vectors
  }

  return(vector_c)
}

# TODO test against simple example
# TODO test that columns sum to one (need to have padded zeroes in both)
# TODO document function and document delay distribution matrix format
#' Title
#'
#' @param vector_a
#' @param matrix_b
#' @param vector_first boolean. Delay described in vector is applied before delay described in matrix
#'
#' @return square matrix. Delay disitribution matrix
.convolve_delay_distribution_vector_with_matrix <- function(vector_a, matrix_b, vector_first = TRUE){

   #TODO add check that matrix_b is square
   N <- nrow(matrix_b)
   # Right-pad vector with zeroes to bring to same dimension as square matrix
   vector_a <- c(vector_a, rep(0, times = max(0, N - length(vector_a))))

   # Initialize result matrix
   convolved_matrix <- matrix(0, nrow = N, ncol = N)

   # Iterate over columns (each column represents the delay distribution on a specific date)
   for(j in 1:N) {
     # Iterate over rows
      for(i in 0 : (N - j)) {
        if(vector_first) { # Take corresponding row in matrix_b
          # The row is left-truncated (only j to N indices) so as to start at same date (date with index j) as column in convolved matrix
          matrix_b_elements <- matrix_b[i + j, j : (j + i) ]
        } else { # Take corresponding column in matrix_b (and revert it)
          matrix_b_elements <- matrix_b[(i + j) : j, j]
        }

        truncated_vector_a <- vector_a[1:i]
        convolved_matrix[i + j, j] <- truncated_vector_a %*% row_matrix_b
      }
   }
   return(convolved_matrix)
}

# TODO test on simple example
# TODO test that columns sum to one (need to have padded zeroes in both)
# TODO document function and document delay distribution matrix format
#' Title
#'
#' Note that this convolution operation is not commutative!
#' @param matrix_a square numeric matrix
#' @param matrix_b square numeric matrix
#'
#' @return square matrix. Convolved matrix of time-varying delays.
.convove_delay_distribution_matrices <- function(matrix_a, matrix_b){
  #TODO return error if matrices are not square or not of the same size

  N <- nrow(matrix_a)
  # Initialize result matrix
  convolved_matrix <- matrix(0, nrow = N, ncol = N)

  # Iterate over columns (each column represents the delay distribution on a specific date)
  for(j in 1:N) {
    # Iterate over rows
    for(i in 0 : (N - j)) {

      # Take truncated column of matrix_a (first delay applied)
      matrix_a_elements <- matrix_b[(i + j) : j, j ]
      # Take truncated row of matrix_b (second delay applied)
      matrix_b_elements <- matrix_b[i + j, j : (j + i) ]

      convolved_matrix[i + j, j] <- matrix_a_elements %*% matrix_b_elements
    }
  }
  return(convolved_matrix)

}

#' Build an empirical Cumulative Distribution Function
#' from the convolution of a gamma distribution and an empirical distribution
#'
#' Utility function that sums samples from an empirical distribution and a gamma distribution.
#'
#' @param shape numeric vector
#' @param scale numeric vector
#' @param number_of_samples integer. Number of samples.
#'
#' @return ecdf object
.make_ecdf_from_empirical_data_and_gamma <- function(gamma_draws,
                                                    empirical_distr,
                                                    multiplier_init = 100){

  multiplier <- multiplier_init
  while(length(gamma_draws) < (length(empirical_distr)*multiplier) && multiplier > 1) {
    multiplier <- floor(multiplier * 0.8)
  }

  if(multiplier < 1) {
    multiplier <- 1
  }

  if(multiplier == 1) {
    final_length <- min(length(gamma_draws), length(empirical_distr))
    draws <- sample(gamma_draws, final_length, replace = F) + sample(empirical_distr, final_length, replace = F)
  } else {
    draws <- gamma_draws[1:(length(empirical_distr)*multiplier)] + rep(empirical_distr, times=multiplier)
  }

  return(stats::ecdf(draws))
}

#TODO fill in and update documentation
#' Build a waiting time distribution from the convolution of two gamma distributions
#'
#' @param parm1_incubation numeric. Shape parameter of the gamma distribution representing the delay between infection and symptom onset.
#' @param parm2_incubation numeric. Scale parameter of the gamma distribution representing the delay between infection and symptom onset.
#' @param parm1_onset_to_report numeric. Shape parameter of the gamma distribution representing the delay between symtom onset and observation.
#' @param parm2_onset_to_report numeric. Scale parameter of the gamma distribution representing the delay between symtom onset and observation.
#' @param distribution_type_incubation string.
#' @param distribution_type_onset_to_report string.
#' @param max_quantile numeric. between 0 and 1.
#'
#' @return vector specifying the CDF between each time step of the waiting time distribution.
#' @export
#'
#' @examples
#' #TODO add examples
combine_incubation_with_reporting_delay <- function(parm1_incubation,
                                                    parm2_incubation,
                                                    parm1_onset_to_report,
                                                    parm2_onset_to_report,
                                                    distribution_type_incubation = "gamma",
                                                    distribution_type_onset_to_report = "gamma",
                                                    max_quantile = 0.9999) {


  delay_distribution_incubation <- build_delay_distribution(parm1 = parm1_incubation,
                                                            parm2 = parm2_incubation,
                                                            distribution_type = distribution_type_incubation,
                                                            max_quantile = max_quantile)

  delay_distribution_onset_to_report <- build_delay_distribution(parm1 = parm1_onset_to_report,
                                                                 parm2 = parm2_onset_to_report,
                                                                 distribution_type = distribution_type_onset_to_report,
                                                                 max_quantile = max_quantile)


  convolved_output <- .convolve_delay_distribution_vectors(delay_distribution_incubation,
                                                           delay_distribution_onset_to_report)

  return(convolved_output)
}




#TODO improve function documentation
#TODO format of empirical_delays must be specified somewhere:
# use "event_date" and "report_delay" as column names

#' Build matrix of delay distributions through time from empirical delay data.
#'
#' This matrix is required for the application of the Richardson-Lucy algorithm.
#'
#' @param empirical_delays tibble. format to be specified
#' @param start_date Date. First date of incidence data
#' @param N integer. Length of incidence time series
#' @param time_step string. "day", "X days", "week", "month"... (see \link{\code{seq.Date}} for details)
#' @param min_number_cases integer. Minimal number of cases to build empirical distribution from
#' @param upper_quantile_threshold numeric. Between 0 and 1. TODO add details
#'
#' @return
#' @export
#'
#' @examples
#' #TODO add example
get_matrix_empirical_waiting_time_distr <- function(empirical_delays,
                                                    start_date,
                                                    N,
                                                    time_step = "day",
                                                    min_number_cases = 300,
                                                    upper_quantile_threshold = 0.99){

  ##TODO need to account for offset if onset data (or not onset data?)
  ##TODO reconsider if we make gamma fit (allow to turn it off, or to use different distribution)

  all_dates <- seq.Date(from = start_date, by = time_step, length.out = N)

  empirical_delays <- empirical_delays %>%
    dplyr::filter(event_date %in% all_dates)

  delay_counts <- empirical_delays %>%
    dplyr::select(report_delay) %>%
    dplyr::group_by(report_delay) %>%
    dplyr::summarise(counts = n(), .groups = "drop")

  #TODO rethink way of defining this right_truncation threshold
  threshold_right_truncation <- delay_counts %>%
    dplyr::mutate(cumul_freq = cumsum(counts)/sum(counts)) %>%
    dplyr::filter(cumul_freq > upper_quantile_threshold) %>%
    utils::head(n=1) %>%
    dplyr::pull(report_delay)

  min_number_cases <- min(min_number_cases, sum(delay_counts$counts))

  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)

  # Populate the delay_distribution_matrix by column
  for(i in 1:N) {

    if(i > (N - threshold_right_truncation) & N > threshold_right_truncation ) {

      delay_distribution_matrix[, i ] <-  c(0, delay_distribution_matrix[1:(N-1), i - 1 ])
      next
    }

    #TODO reconsider if we do week averaging: maybe not
    #TODO get the last X delays (or get the last from the same date)
    weeks_averaged <- 0
    repeat{
      weeks_averaged <- weeks_averaged + 1
      recent_counts_distribution <- empirical_delays %>%
        dplyr::filter( event_date %in% .get_dates_to_average_over(i, all_dates, weeks_averaged)) #TODO: add .get_dates_to_average_over

      if(nrow( recent_counts_distribution ) >= min_number_cases) {
        break
      }
    }

    recent_delays <- recent_counts_distribution %>% dplyr::pull(delay)

    gamma_fit <- try(fitdistrplus::fitdist(recent_delays + 1, distr = "gamma"))
    if ("try-error" %in% class(gamma_fit)) {
      #TODO only output this if verbose output
      cat("    mle failed to estimate the parameters. Trying method = \"mme\"\n")
      gamma_fit <- fitdistrplus::fitdist(recent_delays + 1, distr = "gamma", method = "mme")
    }
    #TODO if none work revert to empirical distribution

    shape_fit <- gamma_fit$estimate["shape"]
    rate_fit <- gamma_fit$estimate["rate"]


    last_index <- N - i + 1
    x <- (1:last_index) + 0.5
    x <- c(0, x)

    cdf_values <- stats::pgamma(x, shape = shape_fit, rate = rate_fit)
    freq <- diff(cdf_values)

    if(length(freq) >= last_index) {
      delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), freq[1:last_index])
    } else {
      delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), freq[1:length(freq)], rep(0, times = last_index - length(freq)))
    }
  }

  return( delay_distribution_matrix )
}

#####

#TODO write functions to convolve two delay distributions together and represent them as functions. This should make the deconvolution easier.



#####
#TODO make additional function that prepares incidence if it is onset data to take into account the fact that it needs to first be reported
.build_delay_distribution_matrix_from_empirical_data <- function(empirical_delay_data) {

  ##TODO finish
  # this is unfinished work
  if(! is_onset_data ) { # copy pasted, maybe we don't keep the if statement
    delay_distribution_matrix_onset_to_report <- .get_matrix_empirical_waiting_time_distr(
      empirical_delays,
      all_dates[(days_further_in_the_past_incubation + 1):length(all_dates)])

    delay_distribution_matrix_incubation <- .get_matrix_constant_waiting_time_distr(
      constant_delay_distribution_incubation,
      all_dates)

    #TODO round this
    initial_delta_incubation <- min(which(cumsum(constant_delay_distribution_incubation) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
    initial_delta_report <-  stats::median(empirical_delays$delay, na.rm = T)
  }

  return()
}


