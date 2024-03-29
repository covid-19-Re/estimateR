% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateR.R
\name{module_methods}
\alias{module_methods}
\title{Methods available for each module}
\arguments{
\item{smoothing_method}{string. Method used to smooth the original incidence data.
Available options are:
\itemize{
\item{'LOESS', implemented in \code{\link{.smooth_LOESS}}}
}}

\item{deconvolution_method}{string. Method used to infer timings of infection
events from the original incidence data (aka deconvolution step).
Available options are:
\itemize{
\item{'Richardson-Lucy delay distribution',
implemented in \code{\link{.deconvolve_incidence_Richardson_Lucy}}}
}}

\item{estimation_method}{string. Method used to estimate reproductive number
values through time from the reconstructed infection timings.
Available options are:
\itemize{
\item{'EpiEstim sliding window',
implemented in \code{\link{.estimate_Re_EpiEstim_sliding_window}}}
\item{'EpiEstim piecewise constant',
implemented in \code{\link{.estimate_Re_EpiEstim_piecewise_constant}}}
}}

\item{uncertainty_summary_method}{string. One of the following options:
\itemize{
\item{'NONE' if no summary of bootstrap estimates is required}
\item{'original estimate - CI from bootstrap estimates'.
The confidence interval is built using bootstrapped estimates
and centered around the original estimates.}
\item{'bagged mean - CI from bootstrap estimates'.
The confidence interval is built using bootstrapped estimates
and centered around the mean of bootstrapped estimates and original estimates.}
}}

\item{bootstrapping_method}{string. Method to perform bootstrapping
of the original incidence data.
Available options are:
\itemize{
\item{'non-parametric block boostrap',
implemented in \code{\link{.block_bootstrap}}}
}}
}
\description{
Methods available for each module
}
