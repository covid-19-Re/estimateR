% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-validation.R
\name{.is_valid_computation_ready_delay_object}
\alias{.is_valid_computation_ready_delay_object}
\title{Utility functions for input validity.}
\usage{
.is_valid_computation_ready_delay_object(
  delay_object,
  parameter_name,
  incidence_data_length
)
}
\arguments{
\item{delay_object}{user input object to be tested}

\item{parameter_name}{A string containing the name the tested parameter had in the initial function in which it was passed.}

\item{incidence_data_length}{A number representing the length of the given incidence data.}
}
\value{
TRUE if all tests were passed. Throws an error otherwise.
}
\description{
Utility function that checks whether a user input is a valid computation-ready delay object.
This means it can be one of the following:
\itemize{
\item a probability distribution vector: a numeric vector with no \code{NA} or negative values, whose entries sum up to 1
\item a delay distribution matrix (as described in \code{\link{.check_is_delay_distribution_matrix}})
}
}
