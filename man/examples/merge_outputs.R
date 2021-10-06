shape_incubation = 3.2
scale_incubation = 1.3
delay_incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

smoothed_incidence <- smooth_incidence(HK_incidence_data$onset_incidence)
deconvolved_incidence <- deconvolve_incidence(
  smoothed_incidence,
  delay = delay_incubation
)

## Basic usage of merge_outputs

merged_incidence_1 <- merge_outputs(
  list("smoothed symptom onset" = smoothed_incidence,
       "deconvolved symptom onset" = deconvolved_incidence)
)


## Advanced usage of merge_outputs

merged_incidence_2 <- merge_outputs(
  list("smoothed symptom onset" = smoothed_incidence,
       "deconvolved symptom onset" = deconvolved_incidence),
  ref_date = HK_incidence_data$date[1],
  include_index = TRUE,
  index_col = "index"
)
