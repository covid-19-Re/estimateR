% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-simulate.R
\name{.random_round}
\alias{.random_round}
\title{Round a value to the integer to either the floor or ceiling value, based on a random draw.}
\usage{
.random_round(observations)
}
\arguments{
\item{observations}{A vector of numeric values.}
}
\value{
Randomly-rounded observations: vector of integers.
}
\description{
The probability of rounding to the ceiling value is equal to
the difference between the unrounded value and its floor value.
}
\details{
For instance, 1.3 has 0.3 probability of being rounded to 2,
and 0.7 probability of being rounded to 1.
}
