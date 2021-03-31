README
================

# *estimateR* package

The *estimateR* package provides tools to estimate the effective
reproductive number through time from delayed and indirect observations
of infection events. This package is under development. The README file
is still at a super rough stage.

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


## The numbers below are part of the user input, 
## they need to be adapted to the particular disease studied.

# Incubation period - gamma distribution parameters
shape_incubation = 3.2 
scale_incubation = 1.3

# Delay from onset of symptoms to case observation - gamma distribution parameters
shape_onset_to_report = 2.7
scale_onset_to_report = 1.6

# We specify these parameters in the same unit as the time steps in the original observation data.
# For instance, if the original data represents daily reports,
# the parameters below must be specified in days (this is the case in this toy example).
mean_serial_interval = 4.8
std_serial_interval = 2.3

# The estimation window corresponds to the size of the sliding window used in EpiEstim
# (see help(estimate_Re) for additional details.
# Here, it is set to three days.
estimation_window = 3 

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
#> 1  2020-02-16                 NA                 NA              18.46059
#> 2  2020-02-17                 NA                 NA              22.01742
#> 3  2020-02-18                 NA                 NA              25.82375
#> 4  2020-02-19                 NA                 NA              30.36457
#> 5  2020-02-20                 NA                 NA              36.13353
#> 6  2020-02-21                 NA                 NA              43.49776
#> 7  2020-02-22                 NA                 NA              52.50639
#> 8  2020-02-23                 NA                 NA              63.27568
#> 9  2020-02-24                  4           17.75345              75.83544
#> 10 2020-02-25                  9           21.35431              90.04614
#> 11 2020-02-26                 19           25.29413             105.45573
#> 12 2020-02-27                 14           30.08292             121.81231
#> 13 2020-02-28                 36           36.23071             139.23459
#> 14 2020-02-29                 16           44.08374             157.86473
#> 15 2020-03-01                 39           53.55827             177.84684
#> 16 2020-03-02                 27           64.51933             196.20617
#> 17 2020-03-03                 46           76.83199             210.77829
#> 18 2020-03-04                 77           90.36127             222.95087
#> 19 2020-03-05                 78          104.73262             234.19047
#> 20 2020-03-06                113          119.79569             246.06494
#>      R_mean
#> 1        NA
#> 2        NA
#> 3        NA
#> 4        NA
#> 5  4.453265
#> 6  3.173586
#> 7  2.662685
#> 8  2.441576
#> 9  2.336766
#> 10 2.273414
#> 11 2.214492
#> 12 2.145094
#> 13 2.065283
#> 14 1.982720
#> 15 1.905197
#> 16 1.826182
#> 17 1.735341
#> 18 1.631936
#> 19 1.528346
#> 20 1.440684
```
