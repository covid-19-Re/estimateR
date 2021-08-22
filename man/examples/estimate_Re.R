## Building incidence_data
# estimate_Re assumes incidence_data represents infections, not delayed noisy
# observations of infections. Thus, we need to first smooth the incidence data
# and then perform a deconvolution step. For more details, see the smooth_incidence
# and deconvolve_incidence functions.

shape_incubation <- 3.2
scale_incubation <- 1.3
delay_incubation <- list(name = "gamma",
                         shape = shape_incubation,
                         scale = scale_incubation)
#todo -> i don t think you use case_incidence here; 
smoothed_incidence <- smooth_incidence(EST_incidence_data$case_incidence)
deconvolved_incidence <- deconvolve_incidence(smoothed_incidence, 
                                              delay_incubation = delay_incubation)


## Basic usage of estimate_Re

Re_estimate <- estimate_Re(incidence_data = deconvolved_incidence)
