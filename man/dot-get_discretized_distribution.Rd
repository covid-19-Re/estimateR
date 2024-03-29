% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-distribution.R
\name{.get_discretized_distribution}
\alias{.get_discretized_distribution}
\title{Discretize a probability distribution.}
\usage{
.get_discretized_distribution(distribution, right_boundary, offset_by_one)
}
\arguments{
\item{distribution}{list. probability distribution specified in list format
e.g. list(name = "gamma", shape = 2, scale = 4).
The \code{distribution} list must contain a 'name' element, this element must be  a string and
correspond to one of the types of \code{\link[stats:Distributions]{distributions}} supported in the \link[stats]{Distributions} package.
\code{distribution} must also contain parameters for the specified distribution, in the form '\code{parameter_name=parameter_value}'.}

\item{right_boundary}{positive numeric value.
Maximum number of time steps to discretize the \code{distribution} over.}

\item{offset_by_one}{boolean.
Set to TRUE if \code{distribution} represents the fit of data that was offset by one
(\code{fitted_data = original_data + 1}) to accommodate zeroes in \code{original_data}.}
}
\value{
vector containing weights of the discretized probability distribution.
}
\description{
Discretize a probability distribution.
}
