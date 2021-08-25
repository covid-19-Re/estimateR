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
delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

total_delay <- convolve_delays(list(delay_incubation, delay_onset_to_report))

deconvolved_incidence <- deconvolve_incidence(
  incidence_data = smoothed_incidence,
  delay = total_delay
)