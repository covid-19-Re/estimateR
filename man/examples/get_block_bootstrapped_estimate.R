## Basic usage of get_block_bootstrapped_estimate
# (Only 10 bootstrap replicates are generated to keep the code fast. In practice,
# use more.)

shape_incubation = 3.2
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma",
                              shape = shape_onset_to_report,
                              scale = scale_onset_to_report)


Re_estimate_1 <- get_block_bootstrapped_estimate(
  HK_incidence_data$case_incidence,
  N_bootstrap_replicates = 10,
  delay = list(delay_incubation, delay_onset_to_report)
)


## Advanced usage of get_block_bootstrapped_estimate
# (Only 10 bootstrap replicates are generated to keep the code fast. In practice,
# use more.)


# Incorporating prior knowledge over Re. Here, Re is assumed constant over a time
# frame of one week, with a prior mean of 1.25.

Re_estimate_2 <- get_block_bootstrapped_estimate(
  HK_incidence_data$case_incidence,
  N_bootstrap_replicates = 10,
  delay = list(delay_incubation, HK_delay_data),
  ref_date = HK_incidence_data$date[1],
  estimation_method = 'EpiEstim piecewise constant',
  interval_length = 7,
  uncertainty_summary_method = 'bagged mean - CI from bootstrap estimates',
  mean_Re_prior = 1.25
)

# Incorporating prior knowledge over the disease. Here, we assume the mean of the
# serial interval to be 5 days, and the deviation is assumed to be 2.5 days. The
# delay between symptom onset and case confirmation is passed as empirical data.

Re_estimate_3 <- get_block_bootstrapped_estimate(
  HK_incidence_data$case_incidence,
  N_bootstrap_replicates = 10,
  delay = list(delay_incubation, HK_delay_data),
  ref_date = HK_incidence_data$date[1],
  mean_serial_interval = 5,
  std_serial_interval = 2.5
)
