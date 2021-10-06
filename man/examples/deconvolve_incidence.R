smoothed_onset_incidence <- smooth_incidence(HK_incidence_data$onset_incidence)
smoothed_case_incidence <- smooth_incidence(HK_incidence_data$case_incidence)

## Deconvolving symptom onset data.
# In case the data to be deconvolved represents noisy observations of symptom
# onset, only the delay distribution of the incubation time needs to be specified
# (time that passes between case incidence and showing of symptoms).

shape_incubation = 3.2
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

deconvolved_incidence_1 <- deconvolve_incidence(
  incidence_data = smoothed_onset_incidence,
  delay = delay_incubation
)


## Deconvolving report incidence data.
# In case the data to be deconvolved represents noisy observations of case reports,
# both the delay distribution of the incubation time and the delay distribution
# of the time that passes between symptom onset and the case being reported.

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma",
                              shape = shape_onset_to_report,
                              scale = scale_onset_to_report)

deconvolved_incidence_2 <- deconvolve_incidence(
  incidence_data = smoothed_case_incidence,
  delay = list(delay_incubation, delay_onset_to_report)
)


## Other available formats for specifying delay distributions

# Discretized delay distribution vector
mean_incubation = 5.2
std_incubation = 1.6
delay_distribution_incubation <- list(name="norm",
                                      mean = mean_incubation,
                                      sd = std_incubation)
delay_incubation_vector <- build_delay_distribution(delay_distribution_incubation)

deconvolved_incidence_3 <- deconvolve_incidence(
  incidence_data = smoothed_onset_incidence,
  delay = delay_incubation_vector
)

# Discretized delay distribution matrix
delay_distribution_matrix <- get_matrix_from_empirical_delay_distr(
  HK_delay_data,
  n_report_time_steps = length(smoothed_case_incidence)
)
deconvolved_incidence_4 <- deconvolve_incidence(
  incidence_data = smoothed_case_incidence,
  delay = list(delay_incubation, delay_distribution_matrix)
)

# Dataframe containing empirical delay data
deconvolved_incidence_5 <- deconvolve_incidence(
  incidence_data = smoothed_case_incidence,
  delay = list(delay_incubation, HK_delay_data)
)


