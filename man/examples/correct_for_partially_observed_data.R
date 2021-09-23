## Basic usage of correct_for_partially_observed_data

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

corrected_incidence_data_1 <- correct_for_partially_observed_data(
  incidence_data = HK_incidence_data$onset_incidence,
  delay_until_final_report = delay_onset_to_report
)


## Advanced usage of correct_for_partially_observed_data
# Only taking into account cases that have a chance of being observed greater
# than 25%. Here, the delay between symptom onset and report is given as 
# empirical delay data, hence it is needed to specify the date of the first
# entry in incidence_data

corrected_incidence_data_2 <- correct_for_partially_observed_data(
  incidence_data = HK_incidence_data$onset_incidence,
  delay_until_final_report = HK_delay_data,
  ref_date = HK_incidence_data$date[1],
  cutoff_observation_probability = 0.25
)