shape_incubation = 3.2 
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)


## Basic usage of estimate_from_combined_observations
Re_estimate_1 <- estimate_from_combined_observations(
  partially_delayed_incidence = HK_incidence_data$onset_incidence,
  fully_delayed_incidence = HK_incidence_data$report_incidence,
  partial_observation_requires_full_observation = TRUE,
  delay_until_partial = delay_incubation,
  delay_until_final_report = delay_onset_to_report
)


## Advanced usage of estimate_from_combined_observations

# Getting a more verbose result. Adding a date column and returning intermediate
# results as well as the Re estimate.
Re_estimate_2 <- estimate_from_combined_observations(
  partially_delayed_incidence = HK_incidence_data$onset_incidence,
  fully_delayed_incidence = HK_incidence_data$report_incidence,
  partial_observation_requires_full_observation = TRUE,
  delay_until_partial = delay_incubation,
  delay_until_final_report = delay_onset_to_report,
  ref_date = HK_incidence_data$date[1],
  output_Re_only = FALSE
)

# Incorporating prior knowledge over Re. Here, Re is assumed constant over a time
# frame of one week, with a prior mean of 1.25.
Re_estimate_3 <- estimate_from_combined_observations(
  partially_delayed_incidence = HK_incidence_data$onset_incidence,
  fully_delayed_incidence = HK_incidence_data$report_incidence,
  partial_observation_requires_full_observation = TRUE,
  delay_until_partial = delay_incubation,
  delay_until_final_report = delay_onset_to_report,
  estimation_method = 'EpiEstim piecewise constant',
  interval_length = 7,
  mean_Re_prior = 1.25
)

# Incorporating prior knowledge over the disease. Here, the mean of the serial 
# interval is assumed to be 5 days, and the standard deviation is assumed to be 
# 2.5 days. 
Re_estimate_4 <- estimate_from_combined_observations(
  partially_delayed_incidence = HK_incidence_data$onset_incidence,
  fully_delayed_incidence = HK_incidence_data$report_incidence,
  partial_observation_requires_full_observation = TRUE,
  delay_until_partial = delay_incubation,
  delay_until_final_report = delay_onset_to_report,
  mean_serial_interval = 5,
  std_serial_interval = 2.5
)