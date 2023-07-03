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
the **[*estimateR*
vignette](https://covid-19-re.github.io/estimateR/articles/estimateR.html)**
and additional vignettes for more details.

We start with arbitrary incidence data representing daily counts of
confirmed cases for a disease of interest. This incidence data counts
the number of people tested positive for a particular infectious disease
and reported on day *t*.

``` r
library(estimateR)
date_first_data_point <- as.Date("2020-02-24")
toy_incidence_data <- c(4,9,19,14,36,16,39,27,46,
                          77,78,113,102,134,165,183,
                          219,247,266,308,304,324,346,
                          348,331,311,267,288,254,239)
```

Below we explain the required user input. It contains assumptions that
need to be adapted to the particular disease and setting studied.

#### Reporting delay

In this toy example, we assume that the incidence data represents cases
by date of report. To infer the time series of cases by date of symptom
onset, we need an assumption about the reporting delay, i.e. the time
between symptom onset and reporting. This delay is specific to each
reporting system and can be estimated from data by public health
authorities (e.g. from line list data). The delay can vary for each case
and is therefore described by a delay distribution.

Here we assume that the delay is Gamma distributed with the following
parameters:

``` r
shape_onset_to_report = 2.7
scale_onset_to_report = 1.6
onset_to_report <- list(name="gamma", shape = shape_onset_to_report, scale = scale_onset_to_report)
```

#### Incubation period

With the reporting delay, we can estimate the time series of cases by
date of symptom onset. However, this is still delayed from the actual
time series of infections. We also need to account for an incubation
period, i.e. the time between infection and symptom onset in patients.
Such information is typically obtained from clinical studies. The
incubation period varies between patients and is therefore described by
a distribution.

Here we assume that the incubation period is Gamma distributed with the
following parameters:

``` r
shape_incubation = 3.2 
scale_incubation = 1.3
incubation <- list(name="gamma", shape = shape_incubation, scale = scale_incubation)
```

#### Generation time / serial interval

To estimate the effective reproductive number from the time series of
infections, we need an assumption about the generation interval,
i.e. the time from infection of a primary case to the infection of a
secondary case. In theory, we are interested in the so-called intrinsic
generation interval distribution, which describes the average
infectiousness of a patient over time and does not depend on the
epidemic phase. In practice however, estimates of the serial interval
(time between symptom onset of the primary case and symptom onsets of
the secondary case) is used as a proxy for the generation interval.
Estimates of the serial interval distribution are typically obtained
from household or contact-tracing studies.

> Please note that the validity of using the serial interval as a proxy
> for the generation interval depends on specific assumptions about
> disease progression and infection, which may be more or less violated
> for a certain disease.

In estimateR, we assume that the generation interval / serial interval
is gamma distributed. We supply the mean and standard deviation of this
gamma distribution:

``` r
mean_serial_interval = 4.8
std_serial_interval = 2.3
```

Note that we specify these parameters in the same time unit as the time
steps of the original observation data. For instance, if the original
data represents daily reports (as in this toy example), then the
parameters must be specified in days.

#### Estimation window

The estimation window corresponds to the size of the smoothing window
used in EpiEstim. See `help(estimate_Re)` for additional details. By default,
it is set to three days. The wider the window, the stronger the
smoothing.

``` r
estimation_window = 3
```

#### Estimating the effective reproductive number

Now we have provided all relevant disease- and setting-specific
parameters. With the function
`estimate_Re_from_noisy_delayed_incidence()`, we can compute the
reproductive number through time from the incidence data and our further
specifications.

``` r
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
```

Let’s look at the most recent Re estimates:

``` r
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

- The column `observed_incidence` shows the originally reported case
  counts.
- The column `smoothed_incidence` shows the case counts after smoothing.
- The column `deconvolved_incidence` shows the estimated number of
  infections.
- The column `Re_estimate` shows point estimates for the reproductive
  number.

**Why are there missing values (NAs)?**

We see that no estimates are available for the 8 most recent days - the
earliest available estimate is 9 days delayed. This is because the
reporting delay and incubation period prevent us from inferring the
number of infections for the most recent past. If individuals are
infected today, it will take time until they develop symptoms and are
reported as cases. Thus, cases reported today mostly tell us something
about past infection, and not much about infections today. As a result,
the reproductive number can only be estimated with a certain delay. This
delay may be partially reduced depending on what data is available, and
through the use of nowcasting. See **[this
vignette](https://covid-19-re.github.io/estimateR/articles/comparison_levels_input_details.html)**
for more details.

## How to cite?

When using this package for a publication, please cite:
- our preprint, Scire et al., 2022, medRxiv (medRxiv https://www.medrxiv.org/content/10.1101/2022.06.30.22277095v1)
(This citation will be updated soon with a peer-reviewed publication.)
- Cori et al., 2013, AJE (https://doi.org/10.1093/aje/kwt133). estimateR builds on EpiEstim, the software presented in this publication.
