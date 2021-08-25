README
================

The **estimateR** package provides tools to estimate the effective
reproductive number through time from delayed and indirect observations
of infection events. This package is currently in a beta version.

## Installation

First, make sure you have installed the **devtools** package locally. If
not, run:

    install.packages("devtools")

in RStudio or other.

Then run:

    library(devtools)
    install_github("covid-19-Re/estimateR")

## Documentation

The full documentation is available
**[here](https://covid-19-re.github.io/estimateR/)**.

## Quick example

We demonstrate below a basic use of the **estimateR** package. Check out
the *estimateR* vignette and additional vignettes for more details.

We start with arbitrary incidence data representing daily counts of case
confirmations for a disease of interest. This incidence data counts the
number of people getting tested positive for a particular infectious
disease X on day *T*. With `estimate_Re_from_noisy_delayed_incidence()`,
we compute the reproductive number through time from this incidence data
and specifications of the delays between infection events and case
confirmations as well as the serial interval.

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
incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

# Delay from onset of symptoms to case observation - gamma distribution parameters
shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

# We specify these parameters in the same unit as the time steps in the original observation data.
# For instance, if the original data represents daily reports,
# the parameters below must be specified in days (this is the case in this toy example).
mean_serial_interval = 4.8
std_serial_interval = 2.3

# The estimation window corresponds to the size of the sliding window used in EpiEstim.
# See help(estimate_Re) for additional details.
# Here, it is set to three days.
estimation_window = 3 

## End of user input
toy_estimates <- estimate_Re_from_noisy_delayed_incidence(toy_incidence_data,
  smoothing_method = "LOESS",
  deconvolution_method = "Richardson-Lucy delay distribution",
  estimation_method = "EpiEstim sliding window",
  delay = list(incubation, onset_to_report),
  estimation_window = estimation_window,
  mean_serial_interval = mean_serial_interval,
  std_serial_interval  = std_serial_interval,
  output_Re_only = FALSE,
  ref_date = date_first_data_point,
  time_step = "day"
)

tail(toy_estimates, n = 20)
#>          date observed_incidence smoothed_incidence deconvolved_incidence
#> 19 2020-03-05                 78           104.7215              234.0993
#> 20 2020-03-06                113           119.7628              245.9896
#> 21 2020-03-07                102           135.7005              257.5132
#> 22 2020-03-08                134           152.7363              267.2654
#> 23 2020-03-09                165           171.0718              276.2164
#> 24 2020-03-10                183           187.9392              285.2313
#> 25 2020-03-11                219           201.4078              295.1153
#> 26 2020-03-12                247           212.9347              305.0930
#> 27 2020-03-13                266           223.9770              314.2166
#> 28 2020-03-14                308           235.9919              323.0184
#> 29 2020-03-15                304           247.8836              331.9951
#> 30 2020-03-16                324           258.1096              341.6153
#> 31 2020-03-17                346           267.4566                    NA
#> 32 2020-03-18                348           276.7116                    NA
#> 33 2020-03-19                331           286.6611                    NA
#> 34 2020-03-20                311           296.5506                    NA
#> 35 2020-03-21                267           305.4301                    NA
#> 36 2020-03-22                288           313.7929                    NA
#> 37 2020-03-23                254           322.1324                    NA
#> 38 2020-03-24                239           330.9420                    NA
#>    Re_estimate
#> 19    1.532188
#> 20    1.443083
#> 21    1.374851
#> 22    1.322608
#> 23    1.280009
#> 24    1.244622
#> 25    1.218098
#> 26    1.199757
#> 27    1.186535
#> 28    1.175290
#> 29    1.165038
#> 30    1.157006
#> 31          NA
#> 32          NA
#> 33          NA
#> 34          NA
#> 35          NA
#> 36          NA
#> 37          NA
#> 38          NA
```
