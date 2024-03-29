% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-distribution.R
\name{.get_delay_matrix_from_delay_distributions}
\alias{.get_delay_matrix_from_delay_distributions}
\title{Build delay distribution matrix from a single delay or a list of delay distributions}
\usage{
.get_delay_matrix_from_delay_distributions(distributions, N = 1, ...)
}
\arguments{
\item{distributions}{single distribution or list of distributions,
each element is either a distribution list or discretized probability distribution vector.}

\item{N}{integer. Dimension of output matrix.
Ignored if a list of distributions is provided.}

\item{...}{
  Arguments passed on to \code{\link[=build_delay_distribution]{build_delay_distribution}}
  \describe{
    \item{\code{max_quantile}}{numeric value between 0 and 1.
Upper quantile reached by the last element in the discretized distribution vector.}
    \item{\code{offset_by_one}}{boolean.
Set to TRUE if \code{distribution} represents the fit of data that was offset by one
(\code{fitted_data = original_data + 1}) to accommodate zeroes in \code{original_data}.}
  }}
}
\value{
delay distribution matrix
}
\description{
Build delay distribution matrix from a single delay or a list of delay distributions
}
