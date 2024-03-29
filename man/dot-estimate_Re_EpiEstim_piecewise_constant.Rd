% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_Re.R
\name{.estimate_Re_EpiEstim_piecewise_constant}
\alias{.estimate_Re_EpiEstim_piecewise_constant}
\title{Estimate Re with EpiEstim in a piecewise-constant fashion}
\usage{
.estimate_Re_EpiEstim_piecewise_constant(
  incidence_input,
  import_incidence_input = NULL,
  minimum_cumul_incidence = 12,
  interval_ends = NULL,
  interval_length = 7,
  mean_serial_interval = 4.8,
  std_serial_interval = 2.3,
  mean_Re_prior = 1,
  output_HPD = FALSE
)
}
\arguments{
\item{incidence_input}{Module input object.
List with two elements:
\enumerate{
\item A numeric vector named \code{values}: the incidence recorded on consecutive time steps.
\item An integer named \code{index_offset}: the offset, counted in number of time steps,
by which the first value in \code{values} is shifted compared to a reference time step
This parameter allows one to keep track of the date of the first value in \code{values}
without needing to carry a \code{date} column around.
A positive offset means \code{values} are delayed in the future compared to the reference values.
A negative offset means the opposite.
}}

\item{import_incidence_input}{NULL or module input object.
List with two elements:
\enumerate{
\item A numeric vector named \code{values}: the incidence recorded on consecutive time steps.
\item An integer named \code{index_offset}: the offset, counted in number of time steps,
by which the first value in \code{values} is shifted compared to a reference time step
This parameter allows one to keep track of the date of the first value in \code{values}
without needing to carry a \code{date} column around.
A positive offset means \code{values} are delayed in the future compared to the reference values.
A negative offset means the opposite.
}
If not NULL, this data represents recorded imported cases.
And then \code{incidence_input} represents only local cases.}

\item{minimum_cumul_incidence}{Numeric value.
Minimum number of cumulated infections before starting the Re estimation.
Default is \code{12} as recommended in Cori et al., 2013.}

\item{interval_ends}{Use with \code{estimation_method = "EpiEstim piecewise constant"}
Integer vector. Optional argument.
If provided, \code{interval_ends} overrides the \code{interval_length} argument.
Each element of \code{interval_ends} specifies the right boundary
of an interval over which Re is assumed to be constant for the calculation.
Values in \code{interval_ends} must be integer values corresponding
with the same numbering of time steps as given by \code{incidence_input}.
In other words, \code{interval_ends} and \code{incidence_input},
use the same time step as the zero-th time step.}

\item{interval_length}{Use with \code{estimation_method = "EpiEstim piecewise constant"}
Positive integer value.
Re is assumed constant over steps of size \code{interval_length}.}

\item{mean_serial_interval}{Numeric positive value. \code{mean_si} for \code{\link[EpiEstim]{estimate_R}}}

\item{std_serial_interval}{Numeric positive value. \code{std_si} for \code{\link[EpiEstim]{estimate_R}}}

\item{mean_Re_prior}{Numeric positive value. \code{mean prior} for \code{\link[EpiEstim]{estimate_R}}}

\item{output_HPD}{Boolean. If TRUE, return the highest posterior density interval with the output.}
}
\value{
If \code{output_HPD = FALSE},
value is a module object (a list of the same kind as \code{incidence_input}).
The \code{values} element of the list then contains the Re estimates.
If \code{output_HPD = TRUE}, a list of three module objects is returned.
\itemize{
\item \code{Re_estimate} contains the Re estimates.
\item \code{Re_highHPD} and \code{Re_lowHPD} contain
the higher and lower boundaries of the HPD interval,
as computed by \code{\link[EpiEstim]{estimate_R}}
}
}
\description{
This function returns piecewise-constant Re estimates.
}
