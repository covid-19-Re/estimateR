#### Make square delay distribution matrix from vector of delay distributions
get_matrix_from_constant_waiting_time_distr <- function(waiting_time_distr, N) {

  if(N >= length(waiting_time_distr)) {
    waiting_time_distr <- c(waiting_time_distr, rep(0, times = N - length(waiting_time_distr)))
  }

  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), waiting_time_distr[1:(N - i + 1)])
  }

  return(delay_distribution_matrix)
}

#### Build empirical CDF from draws summing samples from two gamma distributions
make_ecdf_from_two_gammas <- function(shape, scale, numberOfSamples = 1E6) {
  draws <-
    stats::rgamma(numberOfSamples, shape = shape[1], scale = scale[1]) +
    stats::rgamma(numberOfSamples, shape = shape[2], scale = scale[2])
  return(Vectorize(stats::ecdf(draws)))
}

#### Build empirical CDF from draws from a gamma distribution summmed with an empirical distribution
make_ecdf_from_empirical_data_and_gamma <- function(gamma_draws,
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

get_vector_constant_waiting_time_distr <- function(shape_incubation,
                                                   scale_incubation,
                                                   shape_onset_to_report,
                                                   scale_onset_to_report,
                                                   length_out = 200,
                                                   n_random_samples = 1E6) {

  F_h <- make_ecdf_from_two_gammas(shape = c(shape_incubation, shape_onset_to_report),
                                   scale = c(scale_incubation, scale_onset_to_report),
                                   numberOfSamples = n_random_samples)

  f <- Vectorize(function(x){
    if(x < 0) {
      return(0)
    } else if(x < 0.5) {
      return(F_h(0.5))
    } else {
      return(F_h(round(x + 1E-8) + 0.5) - F_h(round(x + 1E-8) - 0.5))
    }
  })

  x <- 0:(length_out - 1)

  return(f(x))
}


#TODO format of empirical_delays must be specified somewhere:
# use "event_date" and "report_delay" as column names
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
    dplyr::pull(delay)

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
      recent_counts_distribution <- onset_to_report_empirical_delays %>%
        dplyr::filter( onset_date %in% get_dates_to_average_over(i, all_dates, weeks_averaged))

      if(nrow( recent_counts_distribution ) >= min_number_cases) {
        break
      }
    }

    recent_delay_counts <-  recent_counts_distribution %>%
      dplyr::select(delay) %>%
      dplyr::group_by(delay) %>%
      dplyr::summarise(counts = n(), .groups = "drop") %>%
      tidyr::complete(delay  = seq(min(delay), max(delay)),
                      fill = list(counts = 0))

    recent_delays <- recent_counts_distribution %>% dplyr::pull(delay)

    gamma_fit <- try(fitdistrplus::fitdist(recent_delays + 1, distr = "gamma"))
    if ("try-error" %in% class(gamma_fit)) {
      cat("    mle failed to estimate the parameters. Trying method = \"mme\"\n")
      gamma_fit <- fitdistrplus::fitdist(recent_delays + 1, distr = "gamma", method = "mme")
    }

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

#####
#TODO make additional function that prepares incidence if it onset data to take into account the fact that it needs to first be reported
build_delay_distribution_matrix_from_empirical_data <- function(empirical_delay_data) {


  ##TODO finish
  if(! is_onset_data ) { # copy pasted, maybe we don't keep the if statement
    delay_distribution_matrix_onset_to_report <- get_matrix_empirical_waiting_time_distr(
      empirical_delays,
      all_dates[(days_further_in_the_past_incubation + 1):length(all_dates)])

    delay_distribution_matrix_incubation <- get_matrix_constant_waiting_time_distr(
      constant_delay_distribution_incubation,
      all_dates)

    #TODO round this
    initial_delta_incubation <- min(which(cumsum(constant_delay_distribution_incubation) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
    initial_delta_report <-  stats::median(empirical_delays$delay, na.rm = T)
  }

  return()
}

##TODO add a function to deal with empirical delay distribution and build delay distribution matrix
## make sure it can deal with Spanish data specificity
