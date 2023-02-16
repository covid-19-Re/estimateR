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


## Below we explain the required user input. It contains assumptions
## that need to be adapted to the particular disease and setting studied.

# Reporting delay:
#
# In this toy example, we assume that the incidence data represents cases by date of report.
# To infer the time series of cases by date of symptom onset, we need an assumption about 
# the reporting delay, i.e. the time between symptom onset and reporting.
# This delay is specific to each reporting system and can be estimated from data by
# public health authorities (e.g. from line list data). The delay can vary for each case
# and is therefore described by a delay distribution.
#
# Here we assume that the delay is Gamma distributed with the following parameters:
shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)

# Incubation period:
# 
# With the reporting delay, we can estimate the time series of cases by date of symptom onset.
# However, this is still delayed from the actual time series of infections. We also need to account for
# an incubation period, i.e. the time between infection and symptom onset in patients.
# Such information is typically obtained from clinical studies. The incubation period varies between
# patients and is therefore described by a distribution.
#
# Here we assume that the incubation period is Gamma distributed with the following parameters:
shape_incubation = 3.2 
scale_incubation = 1.3
incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)

# Generation time / serial interval:
#
# To estimate the effective reproduction number from the time series of infections, we need an assumption
# about the generation interval, i.e. the time from infection of a primary case to the infection of a secondary
# case. In theory, we are interested in the so-called intrinsic generation interval distribution,
# which describes the average infectiousness of a patient over time and does not depend on the epidemic phase.
# In practice however, estimates of the serial interval (time between symptom onset of the primary case and 
# symptom onsets of the secondary case) is used as a proxy for the generation interval. Please note that the
# validity of this approximation depends on specific assumptions about disease progression and infection, which
# may be more or less violated for a certain disease. Estimates of the serial interval distribution are
# typically obtained from household or contact-tracing studies.
#
# In estimateR, we assume that the generation interval / serial interval is gamma distributed and supply
# an assumption about the mean and standard deviation of this gamma distribution.
# Note that we specify these parameters in the same time unit as the time steps of the original observation data.
# For instance, if the original data represents daily reports (as in this toy example), then the parameters below 
# must be specified in days:
mean_serial_interval = 4.8
std_serial_interval = 2.3

# Estimation window
#
# The estimation window corresponds to the size of the smoothing window used in EpiEstim.
# See help(estimate_Re) for additional details.
# Here, it is set to three days. The wider the window, the stronger the smoothing.
estimation_window = 3 

## End of user input

# Now we can estimate Re from the incidence data
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

# Let's look at the most recent Re estimates.
#
# The column 'observed_incidence' shows the originally reported case counts.
# The column 'smoothed_incidence' shows the case counts after smoothing.
# The column 'deconvolved_incidence' shows the estimated number of infections.
# The column 'Re_estimate' shows point estimates for the reproduction number.
#
# We see that no estimates are available for the 8 most recent days - the
# earliest available estimate is 9 days delayed. This is because the reporting
# delay and incubation period prevent us from inferring the number of infections
# for the most recent past. If a patient is infected today, it will take time
# until the individual develops symptoms and is reported as a case. Thus, the
# reproduction number can only be estimated with a certain delay. This can be
# further improved using nowcasting, which is however not a focus of this package.
tail(toy_estimates, n = 20)
#>          date observed_incidence smoothed_incidence deconvolved_incidence
#> 19 2020-03-05                 78           104.7215              233.9893
#> 20 2020-03-06                113           119.7628              245.8158
#> 21 2020-03-07                102           135.7005              257.2442
#> 22 2020-03-08                134           152.7363              266.8594
#> 23 2020-03-09                165           171.0718              275.6173
#> 24 2020-03-10                183           187.9392              284.3650
#> 25 2020-03-11                219           201.4078              293.8851
#> 26 2020-03-12                247           212.9347              303.3842
#> 27 2020-03-13                266           223.9770              311.9033
#> 28 2020-03-14                308           235.9919              319.9628
#> 29 2020-03-15                304           247.8836              328.0520
#> 30 2020-03-16                324           258.1096              336.6381
#> 31 2020-03-17                346           267.4566                    NA
#> 32 2020-03-18                348           276.7116                    NA
#> 33 2020-03-19                331           286.6611                    NA
#> 34 2020-03-20                311           296.5506                    NA
#> 35 2020-03-21                267           305.4301                    NA
#> 36 2020-03-22                288           313.7929                    NA
#> 37 2020-03-23                254           322.1324                    NA
#> 38 2020-03-24                239           330.9420                    NA
#>    Re_estimate
#> 19    1.531783
#> 20    1.442508
#> 21    1.374039
#> 22    1.321471
#> 23    1.278438
#> 24    1.242481
#> 25    1.215222
#> 26    1.195959
#> 27    1.181623
#> 28    1.169092
#> 29    1.157416
#> 30    1.147870
#> 31          NA
#> 32          NA
#> 33          NA
#> 34          NA
#> 35          NA
#> 36          NA
#> 37          NA
#> 38          NA
```
