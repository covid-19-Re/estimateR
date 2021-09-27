## Basic usage of get_matrix_from_empirical_delay_distr
# Obtaining the deconvolved incidence for the full HK incidence data provided in
# the package: obtaining the delay matrix and then using it to recover the 
# deconvolved incidence.

smoothed_incidence <- smooth_incidence(HK_incidence_data$case_incidence)

shape_incubation = 3.2 
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

delay_matrix_1 <- get_matrix_from_empirical_delay_distr(
  HK_delay_data, 
  n_report_time_steps = length(smoothed_incidence)
)

deconvolved_incidence_1 <- deconvolve_incidence(
  incidence_data = smoothed_incidence,
  delay = list(delay_incubation, delay_matrix_1)
)


## Advanced usage of get_matrix_from_empirical_delay_distr
# Obtaining the deconvolved incidence for a section of the HK incidence data 
# provided in the package: computing the delay matrix, fitting gamma distributions 
# to the columns and then using it to recover the deconvolved incidence for the 
# time-frame of interest

smoothed_partial_incidence <- smooth_incidence(HK_incidence_data[30:90,]$case_incidence)

delay_matrix_2 <- get_matrix_from_empirical_delay_distr(
  HK_delay_data, 
  n_report_time_steps = length(smoothed_partial_incidence),
  fit = "gamma",
  ref_date = HK_incidence_data[30,]$date
)

deconvolved_incidence_2 <- deconvolve_incidence(
  incidence_data = smoothed_partial_incidence,
  delay = list(delay_incubation, delay_matrix_2)
)