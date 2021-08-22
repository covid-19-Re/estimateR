##TOASK: am I using $onset_incidence and $report_incidence right? 
smoothed_onset_incidence <- smooth_incidence(HK_incidence_data$onset_incidence)
smoothed_report_incidence <- smooth_incidence(HK_incidence_data$report_incidence)

## Deconvolving symptom onset data. 
# In case the data to be deconvolved represents noisy observations of symptom
# onset, only the delay distribution of the incubation time needs to be specified 
# (time that passes between case incidence and showing of symptoms).

shape_incubation = 3.2#3.2 
scale_incubation = 13 #1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

deconvolved_incidence_1 <- deconvolve_incidence(incidence_data = smoothed_onset_incidence,
                                              delay_incubation = delay_incubation)

## Deconvolving report incidence data.
# In case the data to be deconvolved represents noisy observations of case reports,
# both the delay distribution of the incubation time and the delay distribution 
# of the time that passes between symptom onset and the case being reported.

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

deconvolved_incidence_2 <- deconvolve_incidence(incidence_data = smoothed_report_incidence,
                                              delay_incubation = delay_incubation,
                                              delay_onset_to_report = delay_onset_to_report)

## Available formats for specifying delay distributions
# A list describing a distribution defined in the stats package, as in the examples above.
mean_incubation = 5.2 
std_incubation = 1.6
delay_incubation <- list(name="norm", mean = mean_incubation, sd = std_incubation)

deconvolved_incidence_1 <- deconvolve_incidence(incidence_data = smoothed_onset_incidence,
                                                delay_incubation = delay_incubation)