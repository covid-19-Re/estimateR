## Convolving the delay between infection and onset of symptoms with the delay
# between onset of symptoms and case report to obtain the final delay between
# infection and case report. Using the resulting delay distribution to recover
# the original infection events

smoothed_incidence <- smooth_incidence(HK_incidence_data$case_incidence)

shape_incubation = 3.2
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma",
                              shape = shape_onset_to_report,
                              scale = scale_onset_to_report)

total_delay_1 <- convolve_delays(list(delay_incubation, delay_onset_to_report))

deconvolved_incidence <- deconvolve_incidence(
  incidence_data = smoothed_incidence,
  delay = total_delay_1
)

## Convolving multiple delays
# In this example it is assumed that the delay between infection and case report
# is composed of three delays; the delay between infection and symptom onset,
# the delay between symptom onset and case testing

shape_incubation = 3.2
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

delay_onset_to_test_taken <- list(name = "norm", mean = 5, sd = 2)

delay_test_to_report <- list(name="norm", mean = 2, sd = 0.5)

total_delay_2 <- convolve_delays(list(delay_incubation,
                                      delay_onset_to_test_taken,
                                      delay_test_to_report))

## Convolving delays of multiple types
# Defining the incubation period as a probability vector, and the delay between
# symptom onset and case observation as a delay matrix

delay_incubation <- c(0.01, 0.1, 0.15, 0.18, 0.17, 0.14, 0.11, 0.07, 0.035, 0.020, 0.015)

delay_matrix <- get_matrix_from_empirical_delay_distr(
  HK_delay_data,
  n_report_time_steps = 50
)

total_delay_3 <- convolve_delays(list(delay_incubation, delay_matrix))

