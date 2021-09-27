## Basic usage of make_tibble_from_output

smoothed_incidence <- smooth_incidence(
  incidence_data = HK_incidence_data$case_incidence,
  smoothing_method = "LOESS"
)

smoothed_incidence_tibble_1 <- make_tibble_from_output(smoothed_incidence)


## Advanced usage of make_tibble_from_output

smoothed_incidence_tibble_2 <- make_tibble_from_output(
  output = smoothed_incidence,
  output_name = "incidence",
  ref_date = HK_incidence_data$date[1]
)