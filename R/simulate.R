compute_discretized_infectiousness <- function(k, shapeG=2.73, scaleG=1.39) {
  ### Expression from Cori et al. 2013, Web appendix 11
  w <- k * pgamma(k, shape=shapeG, scale=scaleG) +
    (k-2)* pgamma(k-2, shape=shapeG, scale=scaleG) +
    (-2) * (k-1) * pgamma(k-1, shape=shapeG, scale=scaleG) +
    shapeG * scaleG * (2 * pgamma(k-1, shape=shapeG+1, scale=scaleG) -
                         pgamma(k-2, shape=shapeG+1, scale=scaleG) -
                         pgamma(k, shape=shapeG+1, scale=scaleG))

  return(w)
}

get_infectiousness_profile <- function(mean_SI = 4.8, sd_SI = 2.3, stopping_threshold = 0.999999){

  shapeG <- (mean_SI - 1)^2 / sd_SI^2
  scaleG <- sd_SI^2 / (mean_SI - 1)

  sum_probability_mass <- 0
  day <- 0
  infectiousness_profile <- c()
  while(sum_probability_mass < stopping_threshold) {
    infectiousness_profile[day + 1] <- compute_discretized_infectiousness(day, shapeG, scaleG)
    sum_probability_mass <- sum_probability_mass + infectiousness_profile[day + 1]
    day <- day + 1
  }

  return(infectiousness_profile)
}

draw_It <- function(Rt, incidence, day, infectiousness_profile) {
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
  sampled_infected_incidence <- rpois(1, Rt[day] * summed_infectiousness)
  return(sampled_infected_incidence)
}

#TODO doc
#' Title
#'
#' @param Rt
#' @param imported_infections
#' @param mean_SI
#' @param sd_SI
#'
#' @return
#' @export
simulate_infections <- function(Rt, imported_infections = 1, mean_SI = 4.8, sd_SI = 2.3){

  length_infections <- max(length(Rt), length(imported_infections))

  # Bring vectors to the same length (length_infections) for future additions
  Rt <- c(Rt, rep(0, times = length_infections-length(Rt)))
  imported_infections <- c(imported_infections, rep(0, times = length_infections-length(imported_infections)))

  infectiousness_profile <- get_infectiousness_profile(mean_SI, sd_SI)

  infections <- imported_infections # Initialize with imports

  for(i in 2:length_infections){
    infections[i] = infections[i] + draw_It(Rt = Rt,
                                            incidence = infections,
                                            day = i,
                                            infectiousness_profile = infectiousness_profile)
  }

  return(infections)
}

#TODO doc
#' Title
#'
#' @param infections
#' @param delay
#'
#' @return
#' @export
simulate_delayed_observations <- function(infections, delay){
  total_delay_distribution <- convolve_delays(delays = delay)

  observations <- sapply(1:length(infections), function(x){compute_Ot(infections, total_delay_distribution, day = x)})
  # Add (optional) noise
  observations <- .add_noise(observations,
                                   noise = noise)

  return(observations)
}

# util
compute_Ot <- function(infections, delay_distribution, day) {
  if(day <1) {
    return(0)
  }
  compute_element_in_sum <- function(x) {
    if(x > (length(delay_distribution) - 1) || x >= day || x < 0) {
      return(0)
    } else {
      return(infections[day - x] * delay_distribution[x + 1])
    }
  }
  raw_summed_observations <- sum(sapply(0:(day-1), compute_element_in_sum))

  return(random_round(raw_summed_observations))
}

# util
random_round <- function(observations){
  # We ensure that the returned number of observations is an integer
  # But we don't just round, otherwise we never see values close to zero.
  rounded_observations <- sapply(observations, function(x) {
    floor(x) + stats::rbinom(n = 1, size = 1, prob = x - floor(x)) })

  return(rounded_observations)
}


# Simplified version of noise generation adapted from Huisman et al.
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

  if (noise$type == 'noiseless'){
    observations = observations
  }

  observations = random_round(observations)
  observations[observations < 0] = 0

  return(observations)
}


#TODO doc
#' Title
#'
#' @param infections
#' @param delay_until_partial
#' @param delay_until_final_report
#' @param prob_partial_observation
#'
#' @return
#' @export
simulate_combined_observations <- function(infections, delay_until_partial,
                                           delay_until_final_report,
                                           prob_partial_observation,
                                           noise = list(type = 'noiseless')){
  delay_until_partial <- .get_delay_distribution(delay_until_partial)
  all_partial_observations <- sapply(1:length(infections),
                                     function(x){compute_Ot(infections, delay_until_partial, day = x)})
  sampled_partial_observations <- sapply(all_partial_observations,
                                         function(x){stats::rbinom(n=1,size = x, prob = prob_partial_observation)})

  unsampled_partial_observations <- all_partial_observations - sampled_partial_observations
  unsampled_partial_observations[unsampled_partial_observations < 0] <- 0
  delay_until_final_report <- .get_delay_distribution(delay_until_final_report)
  final_observations <- sapply(1:length(unsampled_partial_observations),
                               function(x){compute_Ot(unsampled_partial_observations, delay_until_final_report, day = x)})

  # Add (optional) noise
  final_observations <- .add_noise(final_observations,
                                   noise = noise)

  .discard_unsampled_full_observations <- function(partial_observations, delay_until_final_report){
    delay_distribution <- .get_delay_distribution(delay_until_final_report)

    sampled_partial_observations <- partial_observations

    cdf <- cumsum(delay_distribution)
    if(cdf[length(cdf)] < 1) {
      cdf <- c(cdf, 1)
    }

    max_negative_index <- length(cdf) - 1

    for (idx in 0:max_negative_index) {
      sampled_partial_observations[length(partial_observations) - idx] <- stats::rbinom(n = 1,
                                                                                        size = partial_observations[length(partial_observations) - idx],
                                                                                        prob = cdf[idx + 1])
    }

    return(sampled_partial_observations)
  }

  sampled_partial_observations <- .discard_unsampled_full_observations(sampled_partial_observations,
                                                                       delay_until_final_report)

  # Add (optional) noise
  sampled_partial_observations <- .add_noise(sampled_partial_observations,
                                             noise = noise)

  return(data.frame(partially_delayed=sampled_partial_observations, fully_delayed = final_observations))
}
