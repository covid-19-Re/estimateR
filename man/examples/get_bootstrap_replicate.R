## Basic usage of get_bootstrap_replicate

bootstrap_replicate_1 <- get_bootstrap_replicate(
  HK_incidence_data$case_incidence
)


## Advanced usage of get_bootstrap_replicate
# Generate a bootstrap replicate of the incidence data, where case numbers are
# allowed to be decimal numbers, and the output is return as a list.

bootstrap_replicate_2 <- get_bootstrap_replicate(
  HK_incidence_data$case_incidence,
  simplify_output = FALSE,
  round_incidence = FALSE
)

