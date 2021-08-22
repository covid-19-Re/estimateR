## Basic usage of smooth_incidence

smoothed_incidence_1 <- smooth_incidence(
  incidence_data = HK_incidence_data$case_incidence,
  smoothing_method = "LOESS"
)


## Advanced usage of smooth_incidence
# Smoothing the incidence using a LOESS window of 15 days, fitting polynomials 
# of degree 2 in the LOESS algorithm, and averaging the initial Re estimate over
# the first 7 days. 

smoothed_incidence_2 <- smooth_incidence(
  incidence_data = HK_incidence_data$case_incidence,
  smoothing_method = "LOESS",
  simplify_output = FALSE,
  data_points_incl = 15, 
  degree = 2,
  initial_Re_estimate_window = 7
)

 