#' Round a value to the integer to either the floor or ceiling value, based on a random draw.
#'
#' The probability of rounding to the ceiling value is equal to
#' the difference between the unrounded value and its floor value.
#'
#' For instance, 1.3 has 0.3 probability of being rounded to 2,
#' and 0.7 probability of being rounded to 1.
#'
#' @param observations A vector of numeric values.
#'
#' @return Randomly-rounded observations: vector of integers.
.random_round <- function(observations){
  # We ensure that the returned number of observations is an integer
  # But we don't just round, otherwise we never see values close to zero.
  rounded_observations <- sapply(observations, function(x) {
    floor(x) + stats::rbinom(n = 1, size = 1, prob = x - floor(x)) })

  return(rounded_observations)
}


#' Compute the discretized infectiousness for a particular time step.
#'
#' The gamma distribution specified is not directly the serial interval but the serial interval
#' minus one, to allow for a gamma distribution (otherwise infectiousness at day 0 could not be handled),
#' see Cori et al. 2013 for more details.
#'
#' @param k Integer. Infex of the time step.
#' @param shapeG shape of the gamma distribution of the serial interval -1.
#' @param scaleG scale of the gamma distribution of the serial interval -1.
#'
#' @return value of the infectiousness profile at time step k.
.compute_discretized_infectiousness <- function(k, shapeG=2.73, scaleG=1.39) {
  ### Expression from Cori et al. 2013, Web appendix 11
  w <- k * stats::pgamma(k, shape=shapeG, scale=scaleG) +
    (k-2)* stats::pgamma(k-2, shape=shapeG, scale=scaleG) +
    (-2) * (k-1) * stats::pgamma(k-1, shape=shapeG, scale=scaleG) +
    shapeG * scaleG * (2 * stats::pgamma(k-1, shape=shapeG+1, scale=scaleG) -
                         stats::pgamma(k-2, shape=shapeG+1, scale=scaleG) -
                         stats::pgamma(k, shape=shapeG+1, scale=scaleG))

  return(w)
}

#' Compute a discretized infectiousness profile from a serial interval distribution.
#'
#' The serial interval is assumed to be gamma-distributed here,
#' following Cori et al. 2013.
#'
#' @inheritParams simulate
#' @param stopping_threshold Numeric value between 0 and 1.
#' Threshold on cumulative sum of infectiousness profile returned.
#' There is little reason to change to something else than the default value.
#'
#' @return Discretized infectiousness profile through time,
#' represented as a numeric vector with the first value being the infectiousness at time step 0
#' and each subsequent value being the infectiousness on the time step after the previous value.
.get_infectiousness_profile <- function(mean_SI = 4.8, sd_SI = 2.3, stopping_threshold = 0.999999){

  shapeG <- (mean_SI - 1)^2 / sd_SI^2
  scaleG <- sd_SI^2 / (mean_SI - 1)

  sum_probability_mass <- 0
  day <- 0
  infectiousness_profile <- c()
  while(sum_probability_mass < stopping_threshold) {
    infectiousness_profile[day + 1] <- .compute_discretized_infectiousness(day, shapeG, scaleG)
    sum_probability_mass <- sum_probability_mass + infectiousness_profile[day + 1]
    day <- day + 1
  }

  return(infectiousness_profile)
}

#' Draw the number of infections for a particular time step.
#'
#' @inheritParams simulate
#' @param incidence Numeric vector. Incidence of infections through time before the current time step.
#' @param day Integer. Index of the current time step.
#' @param infectiousness_profile Numeric vector. Discretized infectiousness values through time.
#'
#' @return Positive integer. Simulated number of infections for the current time step.
.draw_It <- function(Rt, incidence, day, infectiousness_profile) {
  if(length(Rt) < day || day <1) {
    return(0)
  }

  compute_element_in_sum <- function(x) {
    if(x > (length(infectiousness_profile) - 1) || x >= day || x < 1) {
      return(0)
    } else {
      return(incidence[day - x] * infectiousness_profile[x + 1])
    }
  }
  summed_infectiousness <- sum(sapply(1:(day-1), compute_element_in_sum))
  sampled_infected_incidence <- stats::rpois(1, Rt[day] * summed_infectiousness)
  return(sampled_infected_incidence)
}

#' Compute the number of delayed observations of infections at the current time step.
#'
#' @param infections Vector representing infections through time.
#' @param delay_distribution Numeric vector or matrix. Discretized delay distribution represented as a vector (matrix).
#' @param day Integer. index of the current time step.
#'
#' @return Integer. Number of observations made on a particular time step.
.compute_Ot <- function(infections, delay_distribution, day) {
  if(day <1) {
    return(0)
  }

  if(is.matrix(delay_distribution)) {
    if(day <= nrow(delay_distribution)) {
      delay_distribution_vector <- as.vector(delay_distribution[day, day:1])
    } else {
      delay_distribution_vector < c(0)
    }
  } else {
    delay_distribution_vector <- delay_distribution
  }
  compute_element_in_sum <- function(x) {
    if(x > (length(delay_distribution_vector) - 1) || x >= day || x < 0) {
      return(0)
    } else {
      return(infections[day - x] * delay_distribution_vector[x + 1])
    }
  }
  raw_summed_observations <- sum(sapply(0:(day-1), compute_element_in_sum))

  return(.random_round(raw_summed_observations))
}


#' Add noise to a series of observations.
#'
#' @param observations Numeric vector. Series of observations through time.
#' @param noise List specifying the type of noise and its parameters, if applicable.
#'
#' @return Positive integer vector. Noisy observations.
.add_noise <- function(observations, noise = list(type = 'iid_noise_sd', sd = 1)){

  if (noise$type == 'gaussian'){
    mult_noise <- stats::rnorm(length(observations), mean = 1, sd = noise$sd)
    observations = mult_noise * observations
  }

  if (noise$type ==  'iid_noise_sd'){
    mult_noise <- stats::rnorm(length(observations), mean = 0, sd = noise$sd)

    # y  =  mu * residual
    observations = observations * exp(mult_noise)     # so the error is iid log-normal
  }

  if(noise$type == 'autocorrelated'){
    ar_coeffs <- noise$ar
    AR_noise <- stats::arima.sim(model = list(ar = ar_coeffs), n = length(observations), sd = noise$sd)

    observations = observations * exp(as.vector(AR_noise))
  }

  if (noise$type == 'noiseless'){
    # No noise added.
    observations = observations
  }

  observations = .random_round(observations)
  observations[observations < 0] = 0

  return(observations)
}


#' Generate a list of delay distributions with a gradual transition between two input delay distributions
#'
#' The initial and final delay distributions must be parameterized as gamma distributions with shape and scale parameters.
#' The intermediary distributions are parameterized by scales and shapes that are linear interpolations between the
#' initial and final shapes and scales.
#' @param init_delay List. First delay distribution in the output list
#' @param final_delay List. Last delay distribution in hte output list
#' @param n_time_steps Integer. Number of output list
#'
#' @return List of delay distributions (specified as lists)
.build_list_of_gradually_changing_delays <- function(init_delay, final_delay, n_time_steps) {
  #TODO enforce that both delays are distributions specified as lists and are gamma distributions with scale and shape

  init_distrib_name <- init_delay$name
  final_distrib_name <- final_delay$name

  if(init_distrib_name != final_distrib_name | init_distrib_name!= "gamma") {
    stop("init_delay and final_delay must be two gamma distributions.")
  }

  gradual_shapes <- seq(from = init_delay$shape, to = final_delay$shape, length.out = n_time_steps)
  gradual_scales <- seq(from = init_delay$scale, to = final_delay$scale, length.out = n_time_steps)

  list_of_distributions <- lapply(1:n_time_steps, function(i) {
    return(list(name = "gamma", scale = gradual_scales[i], shape = gradual_shapes[i]))
  })

  return(list_of_distributions)
}
