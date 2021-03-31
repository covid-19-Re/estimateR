README
================

# *estimateR* package

This package is under development. The README file is still at a super
rough stage.

## Quick toy example

We demonstrate below a basic use of the *estimateR* package. Check out
the estimateR vignette for more details.

We start with arbitrary incidence data representing daily counts of case
confirmations for a disease of interest. This incidence data counts the
number of people getting tested positive for a particular infectious
disease X on day *T*.

``` r
library(estimateR)
date_first_data_point <- as.Date("2020-02-24")
toy_incidence_data <- c(4,9,19,14,36,16,39,27,46,
                          77,78,113,102,134,165,183,
                          219,247,266,308,304,324,346,
                          348,331,311,267,288,254,239)


## The numbers below are part of the user input, they need to be adapted to the particular disease studied.
# Incubation period - gamma distribution parameters
shape_incubation = 3.2 
scale_incubation = 1.3

# Delay from onset of symptoms to case observation - gamma distribution parameters
shape_onset_to_report = 2.7
scale_onset_to_report = 1.6

# We specify these parameters in the same unit as the time steps in the original observation data. For instance, if the original data represents daily reports, the parameters below must be specified in days (this is the case in this toy example).
mean_serial_interval = 4.8
std_serial_interval = 2.3

estimation_window = 3 # The estimation window corresponds to the size of the sliding window used in EpiEstim (see help(estimate_Re) for additional details. Here, it is set to three days.

## End of user input

# Build a vector containing the delay distribution
delay_distribution <- get_vector_constant_waiting_time_distr(shape_incubation,
                                                            scale_incubation,
                                                            shape_onset_to_report,
                                                            scale_onset_to_report)

toy_estimates <- smooth_deconvolve_estimate(toy_incidence_data,
                                         delay_distribution,
                                         smoothing_method = "LOESS",
                                         deconvolution_method = "Richardson-Lucy delay distribution",
                                         estimation_method = "EpiEstim sliding window",
                                         estimation_window = estimation_window,
                                         mean_serial_interval = mean_serial_interval,
                                         std_serial_interval  = std_serial_interval,
                                         output_Re_only = FALSE,
                                         ref_date = date_first_data_point,
                                         time_step = "day")

head(toy_estimates, n = 20)
#>          date observed_incidence smoothed_incidence deconvolved_incidence
#> 1  2020-02-16                 NA                 NA              18.45610
#> 2  2020-02-17                 NA                 NA              22.00860
#> 3  2020-02-18                 NA                 NA              25.80581
#> 4  2020-02-19                 NA                 NA              30.34882
#> 5  2020-02-20                 NA                 NA              36.11492
#> 6  2020-02-21                 NA                 NA              43.47508
#> 7  2020-02-22                 NA                 NA              52.47845
#> 8  2020-02-23                 NA                 NA              63.24369
#> 9  2020-02-24                  4           17.75345              75.79922
#> 10 2020-02-25                  9           21.35431              90.00460
#> 11 2020-02-26                 19           25.29413             105.40811
#> 12 2020-02-27                 14           30.08292             121.75955
#> 13 2020-02-28                 36           36.23071             139.17814
#> 14 2020-02-29                 16           44.08374             157.80535
#> 15 2020-03-01                 39           53.55827             177.78476
#> 16 2020-03-02                 27           64.51933             196.14027
#> 17 2020-03-03                 46           76.83199             210.70836
#> 18 2020-03-04                 77           90.36127             222.88002
#> 19 2020-03-05                 78          104.73262             234.12165
#> 20 2020-03-06                113          119.79569             245.99790
#>      R_mean
#> 1        NA
#> 2        NA
#> 3        NA
#> 4        NA
#> 5  4.452439
#> 6  3.173346
#> 7  2.662589
#> 8  2.441557
#> 9  2.336807
#> 10 2.273509
#> 11 2.214613
#> 12 2.145223
#> 13 2.065419
#> 14 1.982872
#> 15 1.905365
#> 16 1.826351
#> 17 1.735490
#> 18 1.632058
#> 19 1.528450
#> 20 1.440786
```
