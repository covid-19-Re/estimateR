## Basic usage of get_infections_from_incidence
# Recovering infection events from case incidence data assuming distinct gamma 
# distributions for the delay between infection and symptom onset, and the delay 
# between symptom onset and case reporting.

shape_incubation = 3.2 
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

infections_1 <- get_infections_from_incidence(
  HK_incidence_data$case_incidence,
  delay = list(delay_incubation, delay_onset_to_report)
)


## Advanced usage of get_infections_from_incidence
# Recovering infection events from symptom onset data, assuming the same delay
# distributions as above

infections_2 <- get_infections_from_incidence(
  HK_incidence_data$onset_incidence,
  delay = delay_incubation,
  is_partially_reported_data = TRUE,
  delay_until_final_report = delay_onset_to_report,
  ref_date = HK_incidence_data$date[1]
)
