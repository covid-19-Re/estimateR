% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateR.R
\name{combining_observations}
\alias{combining_observations}
\title{Details on combining observations}
\description{
Details on combining observations
}
\details{
With this function, one can specify two types of delayed observations of
infection events (in the same epidemic). The two incidence records are
passed with the \code{partially_delayed_incidence} and \code{fully_delayed_incidence}.
These two types of delayed observations must not overlap with one another:
a particular infection event should not be recorded in both time series.

If the two sets of observations are completely independent from one another,
meaning that they represents two different ways infection events
can be observed, with two different delays
then set \code{partial_observation_requires_full_observation} to \code{FALSE}.
Note that a particular infection events should NOT be recorded twice:
it cannot be recorded both in \code{partially_delayed_incidence} and in \code{fully_delayed_incidence}.

An alternative use-case is when the two sets of observations are not independent
from one another. For instance, if to record a "partially-delayed" event,
one had to wait to record it as a "fully-delayed" event first.
A typical example of this occurs when recording symptom onset events:
in most cases, you must first wait until a case is confirmed via a positive test result
to learn about the symptom onset event (assuming the case was symptomatic in the first place).
But you typically do not have the date of onset of symptoms
for all cases confirmed (even assumed they were all symptomatic cases).
In such a case, we set the \code{partial_observation_requires_full_observation} flag
to \code{TRUE} and we call the incidence constructed from events of
symptom onset \code{partially_delayed_incidence} and
the incidence constructed from case confirmation events
\code{fully_delayed_incidence}.
The delay from infection to symptom onset events is
specified with the \code{delay_until_partial} argument.
The delay from symptom onset to positive test in this example is
specified with the \code{delay_until_final_report} argument.
Note that, for a particular patient,
if the date of onset of symptom is known, the patient must not be counted again
in the incidence of case confirmation.
Otherwise, the infection event would have been counted twice.
}
