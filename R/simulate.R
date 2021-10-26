#' Simulate a series of infections
#'
#' Perform a simulation of infections through time, based on a reproductive number course,
#' a series of imported cases and a serial interval distribution.
#' 
#' @example man/examples/simulate_infections.R
#'
#' @param imported_infections Positive integer vector.
#' Must be of length at least one. Does not need to be the same length as \code{Rt}.
#' If \code{imported_infections} is of greater length than \code{Rt} then \code{Rt}
#' is padded with zeroes to reach the same length.
#' @inheritParams simulate
#'
#' @return Integer vector. Simulated infections through time.
#' @export
simulate_infections <- function(Rt, imported_infections = 1, mean_SI = 4.8, sd_SI = 2.3){

  length_infections <- max(length(Rt), length(imported_infections))

  # Bring vectors to the same length (length_infections) for future additions
  Rt <- c(Rt, rep(0, times = length_infections-length(Rt)))
  imported_infections <- c(imported_infections, rep(0, times = length_infections-length(imported_infections)))

  infectiousness_profile <- .get_infectiousness_profile(mean_SI, sd_SI)

  infections <- imported_infections # Initialize with imports

  for(i in 2:length_infections){
    infections[i] = infections[i] + .draw_It(Rt = Rt,
                                            incidence = infections,
                                            day = i,
                                            infectiousness_profile = infectiousness_profile)
  }

  return(infections)
}

#' Simulate a series of delayed observations from a series of infections.
#'
#' @example man/examples/simulate_delayed_observations.R
#'
#' @inheritParams simulate
#' @inheritParams delay_high
#'
#' @return Integer vector. Simulated delayed observations.
#' @export
simulate_delayed_observations <- function(infections, delay, noise = list(type = "noiseless")){
  total_delay_distribution <- convolve_delays(delays = delay)

  observations <- sapply(1:length(infections), function(x){.compute_Ot(infections, total_delay_distribution, day = x)})
  # Add (optional) noise
  observations <- .add_noise(observations,
                                   noise = noise)

  return(observations)
}

#' Simulate a series of observations from a course of infections, combining some partially-delayed and some fully-delayed observations.
#'
#' The infections that are observed as partially-delayed observations cannot be observed a second time as fully-delayed observations,
#' meaning that they do not show up a second time in the "fully-delayed" column of the result.
#' However, a partially-delayed observation can only be "registered" (included in the "partially-delayed" column) if
#' it is has been virtually observed as a fully-delayed observation first.
#'
#' @inheritParams simulate
#' @inheritParams delay_high
#' @param prob_partial_observation Numeric value between 0 and 1.
#' Probability of an infection to be observed as a partially-delayed observation,
#' instead of as a fully-delayed observation.
#'
#' @example man/examples/simulate_combined_observations.R
#'
#' @return A dataframe containing two columns:
#' a column "partially_delayed" containing partially-delayed observations
#' and a column "fully_delayed" containing fully-delayed observations.
#' @export
simulate_combined_observations <- function(infections,
                                           delay_until_partial,
                                           delay_until_final_report,
                                           prob_partial_observation,
                                           noise = list(type = 'noiseless')){
  delay_until_partial <- .get_delay_distribution(delay_until_partial)
  all_partial_observations <- sapply(1:length(infections),
                                     function(x){.compute_Ot(infections, delay_until_partial, day = x)})
  sampled_partial_observations <- sapply(all_partial_observations,
                                         function(x){stats::rbinom(n=1,size = x, prob = prob_partial_observation)})

  unsampled_partial_observations <- all_partial_observations - sampled_partial_observations
  unsampled_partial_observations[unsampled_partial_observations < 0] <- 0
  delay_until_final_report <- .get_delay_distribution(delay_until_final_report)
  final_observations <- sapply(1:length(unsampled_partial_observations),
                               function(x){.compute_Ot(unsampled_partial_observations, delay_until_final_report, day = x)})

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
