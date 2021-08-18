test_that("get_bootstrap_replicate outputs difference values with same median and bounds as original difference values", {
  skip_on_cran()
  expect_bootstrapped_diff_bounded_by_original_diff <- function(...) {
    data_points_incl <- 21
    degree <- 1
    original_values <- sample.int(1000, size = 10000, replace = TRUE)

    log_original <- log(original_values + 1)
    smoothed_log <- smooth_incidence(log_original, smoothing_method = "LOESS", data_points_incl = data_points_incl)

    diff_smoothed_original <- log_original - smoothed_log

    bootstrap_replicate <- get_bootstrap_replicate(
      incidence_data = original_values,
      bootstrapping_method = "non-parametric block boostrap",
      data_points_incl = data_points_incl,
      degree = degree,
      round_incidence = FALSE
    )

    log_bootstrap <- log(bootstrap_replicate + 1)
    diff_smoothed_bootstrap <- log_bootstrap - smoothed_log

    expect_equal(median(diff_smoothed_bootstrap), median(diff_smoothed_original), tolerance = 0.1)
    expect_gte(min(diff_smoothed_bootstrap), min(diff_smoothed_original) - 1E-1)
    expect_lte(max(diff_smoothed_bootstrap), max(diff_smoothed_original) + 1E-1)
  }

  sapply(1:10, expect_bootstrapped_diff_bounded_by_original_diff)
})
