## Basic usage of estimate_Re_from_noisy_delayed_incidence
shape_incubation = 3.2 
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

Re_estimate_1 <- estimate_Re_from_noisy_delayed_incidence(
  incidence_data = HK_incidence_data$case_incidence,
  delay = list(delay_incubation, delay_onset_to_report)
)

## Advanced usage of estimate_Re_from_noisy_delayed_incidence
# Incorporating prior knowledge over Re. Here, Re is assumed constant over a time
# frame of one week, with a prior mean of 1.25.
Re_estimate_2 <- estimate_Re_from_noisy_delayed_incidence(
  incidence_data = HK_incidence_data$case_incidence,
  delay = list(delay_incubation, delay_onset_to_report),
  estimation_method = "EpiEstim piecewise constant",
  interval_length = 7,
  mean_Re_prior = 1.25
)

# Incorporating prior knowledge over the disease. Here, the mean of the serial 
# interval is assumed to be 5 days, and the standard deviation is assumed to be 
# 2.5 days. 
Re_estimate_3 <- estimate_Re_from_noisy_delayed_incidence(
  incidence_data = HK_incidence_data$case_incidence,
  delay = list(delay_incubation, delay_onset_to_report),
  mean_serial_interval = 5,
  std_serial_interval = 1.25
)

# Incorporating prior knowledge over the epidemic. Here, it is assumed that Re 
# changes values 4 times during the epidemic, so the intervals over which Re is
# assumed to be constant are passed as a parameter.
last_interval_index <- length(deconvolved_incidence$values) + deconvolved_incidence$index_offset 

Re_estimate_4 <- estimate_Re_from_noisy_delayed_incidence(
  incidence_data = HK_incidence_data$case_incidence,
  delay = list(delay_incubation, delay_onset_to_report),
  interval_ends = c(50, 75, 100, 160, last_interval_index)
)

# Recovering the Re HPD as well.
Re_estimate <- estimate_Re_from_noisy_delayed_incidence(
  incidence_data = HK_incidence_data$case_incidence,
  delay = list(delay_incubation, delay_onset_to_report),
  output_HPD = TRUE
)