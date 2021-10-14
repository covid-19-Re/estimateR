## Basic use of simulate_combined_observations
# Simulating combined observations, assuming two gamma delays between infection
# and symptom onset, and symptom onset and case report respectively. It is assumed
# that 20% of the cases are observed as partially-delayed observations.

Re_evolution <- c(rep(2.3, 100))
incidence <- simulate_infections(Re_evolution)

shape_incubation = 3.2
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma",
                              shape = shape_onset_to_report,
                              scale = scale_onset_to_report)
simulated_combined_observations_1 <- simulate_combined_observations(
  incidence,
  delay_until_partial = delay_incubation,
  delay_until_final_report = delay_onset_to_report,
  prob_partial_observation = 0.2
)

## Advanced use of simulate_combined_observations
# Adding gaussian noise to the combined observations simulated above.
simulated_combined_observations_2 <- simulate_combined_observations(
  incidence,
  delay_until_partial = delay_incubation,
  delay_until_final_report = delay_onset_to_report,
  prob_partial_observation = 0.2,
  noise = list(type = 'gaussian', sd = 0.8)
)

