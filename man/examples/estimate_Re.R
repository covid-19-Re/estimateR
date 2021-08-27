## Building incidence_data
# estimate_Re assumes incidence_data represents infections, not delayed noisy
# observations of infections. Thus, we need to first smooth the incidence data
# and then perform a deconvolution step. For more details, see the smooth_incidence
# and deconvolve_incidence functions.

shape_incubation <- 3.2
scale_incubation <- 1.3
delay_incubation <- list(name = "gamma", shape = shape_incubation, scale = scale_incubation)

shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
delay_onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

smoothed_incidence <- smooth_incidence(HK_incidence_data$case_incidence)
deconvolved_incidence <- deconvolve_incidence(
  smoothed_incidence, 
  delay = list(delay_incubation, delay_onset_to_report)
)


## Basic usage of estimate_Re
Re_estimate_1 <- estimate_Re(incidence_data = deconvolved_incidence)


## Advanced usage of estimate_Re
# Incorporating prior knowledge over Re. Here, Re is assumed constant over a time
# frame of one week, with a prior mean of 1.25.
Re_estimate_2 <- estimate_Re(
  incidence_data = deconvolved_incidence,
  estimation_method = 'EpiEstim piecewise constant',
  interval_length = 7,
  mean_Re_prior = 1.25
)

# Incorporating prior knowledge over the disease. Here, the mean of the serial 
# interval is assumed to be 5 days, and the standard deviation is assumed to be 
# 2.5 days. 
Re_estimate_3 <- estimate_Re(
  incidence_data = deconvolved_incidence,
  mean_serial_interval = 5,
  std_serial_interval = 2.5
)

# Incorporating prior knowledge over the epidemic. Here, it is assumed that Re 
# changes values 4 times during the epidemic, so the intervals over which Re is
# assumed to be constant are passed as a parameter.
last_interval_index <- length(deconvolved_incidence$values) + deconvolved_incidence$index_offset 

Re_estimate_4 <- estimate_Re(
  incidence_data = deconvolved_incidence,
  estimation_method = "EpiEstim piecewise constant",
  interval_ends = c(50, 75, 100, 160, last_interval_index)
)

# Recovering the Re HPD as well.
Re_estimate_5 <- estimate_Re(
  incidence_data = deconvolved_incidence,
  output_HPD = TRUE
)