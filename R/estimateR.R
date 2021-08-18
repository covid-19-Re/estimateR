#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom rlang :=
NULL

#' Utility functions for input validity.
#'
#' @param string_user_input A string containing the value that the user passed for the tested string type parameter.
#' @param parameter_name A string containing the name the tested parameter had in the initial function in which it was passed.
#' @param incidence_data_length A number representing the length of the given incidence data.
#' @param user_inputs A list of lists with two elements: the first is the value of the parameter to be tested. The second is the expected type of that parameter.
#' @param number The value to be tested
#' @param user_input The variable to be tested
#'
#' @return TRUE if all tests were passed. Throws an error otherwise.
#'
#' @name validation_utility_params
NULL

#' Module structure characteristics
#'
#' @param incidence_data An object containing incidence data through time.
#' It can either be:
#' \itemize{
#' \item A list with two elements:
#'  \enumerate{
#'  \item A numeric vector named \code{values}: the incidence recorded on consecutive time steps.
#'  \item An integer named \code{index_offset}: the offset, counted in number of time steps,
#'  by which the first value in \code{values} is shifted compared to a reference time step
#'  This parameter allows one to keep track of the date of the first value in \code{values}
#'  without needing to carry a \code{date} column around.
#'  A positive offset means \code{values} are delayed in the future compared to the reference values.
#'  A negative offset means the opposite.
#'  }
#'  \item A numeric vector. The vector corresponds to the \code{values} element
#'  descrived above, and \code{index_offset} is implicitely zero.
#'  This means that the first value in \code{incidence_data}
#'  is associated with the reference time step (no shift towards the future or past).
#' }
#' @param import_incidence_data NULL or argument with the same requirements as \code{incidence_data}.
#' If not NULL, this argument represents records of imported cases and
#' \code{incidence_data} represents local cases only.
#' @param partially_delayed_incidence An object containing incidence data through time.
#' It can be:
#' \itemize{
#' \item A list with two elements:
#'  \enumerate{
#'  \item A numeric vector named \code{values}: the incidence recorded on consecutive time steps.
#'  \item An integer named \code{index_offset}: the offset, counted in number of time steps,
#'  by which the first value in \code{values} is shifted compared to a reference time step
#'  This parameter allows one to keep track of the date of the first value in \code{values}
#'  without needing to carry a \code{date} column around.
#'  A positive offset means \code{values} are delayed in the future compared to the reference values.
#'  A negative offset means the opposite.
#'  }
#'  \item A numeric vector. The vector corresponds to the \code{values} element
#'  descrived above, and \code{index_offset} is implicitely zero.
#'  This means that the first value in \code{incidence_data}
#'  is associated with the reference time step (no shift towards the future or past).
#' }
#' @param fully_delayed_incidence An object containing incidence data through time.
#' It can be:
#' \itemize{
#' \item A list with two elements:
#'  \enumerate{
#'  \item A numeric vector named \code{values}: the incidence recorded on consecutive time steps.
#'  \item An integer named \code{index_offset}: the offset, counted in number of time steps,
#'  by which the first value in \code{values} is shifted compared to a reference time step
#'  This parameter allows one to keep track of the date of the first value in \code{values}
#'  without needing to carry a \code{date} column around.
#'  A positive offset means \code{values} are delayed in the future compared to the reference values.
#'  A negative offset means the opposite.
#'  }
#'  \item A numeric vector. The vector corresponds to the \code{values} element
#'  descrived above, and \code{index_offset} is implicitely zero.
#'  This means that the first value in \code{incidence_data}
#'  is associated with the reference time step (no shift towards the future or past).
#' }
#' @param simplify_output boolean. Return a numeric vector instead of module
#'  output object if output offset is zero?
#'
#'
#' @return A list with two elements:
#'  \enumerate{
#'  \item A numeric vector named \code{values}: the result of the computations on the input data.
#'  \item An integer named \code{index_offset}: the offset, counted in number of time steps,
#'  by which the result is shifted compared to an \code{index_offset} of \code{0}.
#'  This parameter allows one to keep track of the date of the first value in \code{values}
#'  without needing to carry a \code{date} column around.
#'  A positive offset means \code{values} are delayed in the future compared to the reference values.
#'  A negative offset means the opposite.
#'  Note that the \code{index_offset} of the output of the function call
#'  accounts for the (optional) \code{index_offset} of the input.
#'  }
#'  If \code{index_offset} is \code{0} and \code{simplify_output = TRUE},
#'  the \code{index_offset} is dropped and the \code{values}
#'  element is returned as a numeric vector.
#'
#'
#' @name module_structure
NULL


#' Details on combining observations
#' @details
#' With this function, one can specify two types of delayed observations of
#' infection events (in the same epidemic). The two incidence records are
#' passed with the \code{partially_delayed_incidence} and \code{fully_delayed_incidence}.
#' These two types of delayed observations must not overlap with one another:
#' a particular infection event should not be recorded in both time series.
#'
#' If the two sets of observations are completely independent from one another,
#' meaning that they represents two different ways infection events
#' can be observed, with two different delays
#' then set \code{partial_observation_requires_full_observation} to \code{FALSE}.
#' The \code{delay_until_final_report} delay then corresponds to the delay
#' from infection until observation in \code{fully_delayed_incidence}.
#' And the \code{delay_until_partial} delay then corresponds to the delay
#' from infection until observation in \code{partially_delayed_incidence}.
#' Note that a particular infection events should NOT be recorded twice:
#' it cannot be recorded both in \code{partially_delayed_incidence} and in \code{fully_delayed_incidence}.
#'
#' An alternative use-case is when the two sets of observations are not independent
#' from one another. For instance, if to record a "partially-delayed" event,
#' one had to wait to record it as a "fully-delayed" event first.
#' A typical example of this occurs when recording symptom onset events:
#' in most cases, you must first wait until a case is confirmed via a positive test result
#' to learn about the symptom onset event (assuming the case was symptomatic in the first place).
#' But you typically do not have the date of onset of symptoms
#' for all cases confirmed (even assumed they were all symptomatic cases).
#' In such a case, we set the \code{partial_observation_requires_full_observation} flag
#' to \code{TRUE} and we call the incidence constructed from events of
#' symptom onset \code{partially_delayed_incidence} and
#' the incidence constructed from case confirmation events
#' \code{fully_delayed_incidence}.
#' The full delay from infection to positive test in this example is
#' specified with the \code{delay_until_final_report} argument.
#' The delay from infection to symptom onset events is
#' specified with the \code{delay_until_partial} argument.
#' Note that, for a particular patient,
#' if the date of onset of symptom is known, the patient must not be counted again
#' in the incidence of case confirmation.
#' Otherwise, the infection event would have been counted twice.
#'
#'
#' @name combining_observations
NULL

#' Inner module option characteristics
#'
#' @param incidence_input,input,output,input_a,input_b Module input object.
#' List with two elements:
#'  \enumerate{
#'  \item A numeric vector named \code{values}: the incidence recorded on consecutive time steps.
#'  \item An integer named \code{index_offset}: the offset, counted in number of time steps,
#'  by which the first value in \code{values} is shifted compared to a reference time step
#'  This parameter allows one to keep track of the date of the first value in \code{values}
#'  without needing to carry a \code{date} column around.
#'  A positive offset means \code{values} are delayed in the future compared to the reference values.
#'  A negative offset means the opposite.
#'  }
#'
#' @return Module input object.
#' List with two elements:
#'  \enumerate{
#'  \item A numeric vector named \code{values}: the incidence recorded on consecutive time steps.
#'  \item An integer named \code{index_offset}: the offset, counted in number of time steps,
#'  by which the first value in \code{values} is shifted compared to a reference time step
#'  This parameter allows one to keep track of the date of the first value in \code{values}
#'  without needing to carry a \code{date} column around.
#'  A positive offset means \code{values} are delayed in the future compared to the reference values.
#'  A negative offset means the opposite.
#'  }
#' @name inner_module
NULL

#' EpiEstim wrappers arguments
#'
#' @param minimum_cumul_incidence Numeric value.
#' Minimum number of cumulated infections before starting the Re estimation.
#' Default is \code{12} as recommended in Cori et al., 2013.
#' @param mean_serial_interval Numeric positive value. \code{mean_si} for \code{\link[EpiEstim]{estimate_R}}
#' @param std_serial_interval Numeric positive value. \code{std_si} for \code{\link[EpiEstim]{estimate_R}}
#' @param mean_Re_prior Numeric positive value. \code{mean prior} for \code{\link[EpiEstim]{estimate_R}}
#' @param output_HPD Boolean. If TRUE, return the highest posterior density interval with the output.
#' @param import_incidence_input NULL or module input object.
#' List with two elements:
#'  \enumerate{
#'  \item A numeric vector named \code{values}: the incidence recorded on consecutive time steps.
#'  \item An integer named \code{index_offset}: the offset, counted in number of time steps,
#'  by which the first value in \code{values} is shifted compared to a reference time step
#'  This parameter allows one to keep track of the date of the first value in \code{values}
#'  without needing to carry a \code{date} column around.
#'  A positive offset means \code{values} are delayed in the future compared to the reference values.
#'  A negative offset means the opposite.
#'  }
#'  If not NULL, this data represents recorded imported cases.
#'  And then \code{incidence_input} represents only local cases.
#'
#' @return If \code{output_HPD = FALSE},
#' value is a module object (a list of the same kind as \code{incidence_input}).
#' The \code{values} element of the list then contains the Re estimates.
#' If \code{output_HPD = TRUE}, a list of three module objects is returned.
#' \itemize{
#' \item \code{Re_estimate} contains the Re estimates.
#' \item \code{Re_highHPD} and \code{Re_lowHPD} contain
#' the higher and lower boundaries of the HPD interval,
#' as computed by \code{\link[EpiEstim]{estimate_R}}
#' }
#'
#' @name EpiEstim_wrapper
NULL

#' Universal parameters
#'
#' @param verbose Boolean. Print verbose output?
#'
#' @name universal_params
NULL

#' Pipe parameters
#'
#' @param output_Re_only boolean. Should the output only contain Re estimates?
#' (as opposed to containing results for each intermediate step)
#' @param output_infection_incidence_only boolean.
#' Should the output contain only the estimated infection incidence?
#' (as opposed to containing results for intermediary steps)
#'
#' @name pipe_params
NULL

#' Bootstrapping parameters
#'
#' @param N_bootstrap_replicates integer. Number of bootstrap samples.
#' @param combine_bootstrap_and_estimation_uncertainties boolean.
#' If TRUE, the uncertainty interval reported is the union of
#' the highest posterior density interval from the Re estimation
#' with the confidence interval from the boostrapping of time series
#' of observations.
#'
#' @name bootstrap_params
NULL

#' Bootstrapping pipe
#'
#' @return Effective reproductive estimates through time with confidence interval boundaries.
#' If \code{output_Re_only} is \code{FALSE}, then transformations made
#' on the input observations during calculations are output as well.
#'
#' @name bootstrap_return
NULL

#' High-level delay parameters
#'
#' @param delays List of delays, with flexible structure.
#' Each delay in the \code{delays} list can be one of:
#' \itemize{
#' \item{a list representing a distribution object}
#' \item{a discretized delay distribution vector}
#' \item{a discretized delay distribution matrix}
#' \item{a dataframe containing empirical delay data}
#' }
#' @param delay Single delay or list of delays.
#' Each delay can be one of:
#' \itemize{
#' \item{a list representing a distribution object}
#' \item{a discretized delay distribution vector}
#' \item{a discretized delay distribution matrix}
#' \item{a dataframe containing empirical delay data}
#' }
#' @param delay_distribution_final_report
#' Distribution of the delay between the events collected in the incidence data
#' and the a posteriori observations of these events.
#' @param cutoff_observation_probability value between 0 and 1.
#' Only datapoints for timesteps that have a probability of observing a event
#' higher than \code{cutoff_observation_probability} are kept.
#' The few datapoints with a lower probability to be observed are trimmed off
#' the tail of the timeseries.
#' @param is_partially_reported_data boolean.
#' Set to \code{TRUE} if \code{incidence_data} represents delayed observations
#' of infection events that themselves rely on further-delayed observations.
#' @param delay_until_partial Single delay or list of delays.
#' Each delay can be one of:
#' \itemize{
#' \item{a list representing a distribution object}
#' \item{a discretized delay distribution vector}
#' \item{a discretized delay distribution matrix}
#' \item{a dataframe containing empirical delay data}
#' }
#' @param delay_until_final_report Single delay or list of delays.
#' Each delay can be one of:
#' \itemize{
#' \item{a list representing a distribution object}
#' \item{a discretized delay distribution vector}
#' \item{a discretized delay distribution matrix}
#' \item{a dataframe containing empirical delay data}
#' }
#' @param partial_observation_requires_full_observation boolean
#' Set to \code{TRUE} if \code{partially_delayed_incidence} represent
#' delayed observations of infection events that
#' themselves rely on further-delayed observations.
#' See Details for more details.
#'
#' @name delay_high
NULL

#' Empirical delays parameters
#'
#' @param empirical_delays,delays dataframe containing the empirical data. See Details.
#' @param n_report_time_steps integer. Length of the incidence time series in the accompanying analysis.
#' This argument is needed to determine the dimensions of the output matrix.
#' @param min_number_cases integer.
#' Minimal number of cases to build the empirical distribution from.
#' If \code{num_steps_in_a_unit} is \code{NULL}, for any time step T,
#' the \code{min_number_cases} records prior to T are used.
#' If less than \code{min_number_cases} delays were recorded before T,
#' then T is ignored and the \code{min_number_cases} earliest-recorded delays are used.
#' If \code{num_steps_in_a_unit} is given a value, a similar same procedure is applied,
#' except that, now at least \code{min_number_cases} must be taken over a round number of
#' time units. For example, if \code{num_steps_in_a_unit = 7}, and time steps represent consecutive days,
#' to build the distribution for time step T,
#' we find the smallest number of weeks starting from T and going in the past,
#' for which at least \code{min_number_cases} delays were recorded.
#' We then use all the delays recorded during these weeks.
#' Weeks are not meant as necessarily being Monday to Sunday,
#' but simply 7 days in a row, e.g. it can be Thursday-Wednesday.
#' Again, if less than \code{min_number_cases} delays were recorded before T,
#' then T is ignored.
#' We then find the minimum number of weeks, starting from the first recorded delay
#' that contains at least \code{min_number_cases}.
#'
#' @param min_number_cases_fraction numeric. Between 0 and 1.
#' If \code{min_number_cases} is not provided (kept to \code{NULL}),
#' the number of most-recent cases used to build
#' the instant delay distribution is \code{min_number_cases_fraction}
#' times the total number of reported delays.
#' @param min_min_number_cases numeric. Lower bound
#' for number of cases used to build an instant delay distribution.
#' @param upper_quantile_threshold numeric. Between 0 and 1.
#' Argument for internal use.
#' @param fit string. One of "gamma" or "none". Specifies the type of fit that
#' is applied to the columns of the delay matrix
#' @param date_of_interest Date. Date for which the most recent recorded delays are sought.
#' @param num_steps_in_a_unit Optional argument.
#' Number of time steps in a full time unit (e.g. 7 if looking at weeks).
#' If set, the delays used to build a particular
#' delay distribution will span over a round number of such time units.
#' This option is included for comparison with legacy code.
#'
#'
#' @name delay_empirical
NULL



#' Dating parameters
#'
#' @param ref_date Date. Optional. Date of the first data entry in \code{incidence_data}
#' @param time_step string. Time between two consecutive incidence datapoints.
#' "day", "2 days", "week", "year"... (see \code{\link[base]{seq.Date}} for details)
#'
#' @name dating
NULL

#' Methods available for each module
#'
#' @param smoothing_method string. Method used to smooth the original incidence data.
#' Available options are:
#' \itemize{
#' \item{'LOESS', implemented in \code{\link{.smooth_LOESS}}}
#' }
#' @param deconvolution_method string. Method used to infer timings of infection
#' events from the original incidence data (aka deconvolution step).
#' Available options are:
#' \itemize{
#' \item{'Richardson-Lucy delay distribution',
#' implemented in \code{\link{.deconvolve_incidence_Richardson_Lucy}}}
#' }
#' @param estimation_method string. Method used to estimate reproductive number
#' values through time from the reconstructed infection timings.
#' Available options are:
#' \itemize{
#' \item{'EpiEstim sliding window',
#' implemented in \code{\link{.estimate_Re_EpiEstim_sliding_window}}}
#' \item{'EpiEstim piecewise constant',
#' implemented in \code{\link{.estimate_Re_EpiEstim_piecewise_constant}}}
#' }
#' @param uncertainty_summary_method string. One of the following options:
#' \itemize{
#' \item{'NONE' if no summary of bootstrap estimates is required}
#' \item{'original estimate - CI from bootstrap estimates'.
#' The confidence interval is built using bootstrapped estimates
#' and centered around the original estimates.}
#' \item{'bagged mean - CI from bootstrap estimates'.
#' The confidence interval is built using bootstrapped estimates
#' and centered around the mean of bootstrapped estimates and original estimates.}
#' }
#' @param bootstrapping_method string. Method to perform bootstrapping
#' of the original incidence data.
#' Available options are:
#' \itemize{
#' \item{'non-parametric block boostrap',
#' implemented in \code{\link{.block_bootstrap}}}
#' }
#'
#' @name module_methods
NULL

#' Distribution
#'
#' @param distribution list. probability distribution specified in list format
#' e.g. list(name = "gamma", shape = 2, scale = 4).
#' The \code{distribution} list must contain a 'name' element, this element must be  a string and
#' correspond to one of the types of \code{\link[stats:Distributions]{distributions}} supported in the \link[stats]{Distributions} package.
#' \code{distribution} must also contain parameters for the specified distribution, in the form '\code{parameter_name=parameter_value}'.
#' @param vector_a,vector_b,delay_distribution_vector discretized probability distribution vector
#' @param matrix_a,matrix_b,delay_distribution_matrix discretized delay distribution matrix
#' @param quantile Value between 0 and 1. Quantile of the distribution.
#'
#' @name distribution
NULL

#' Empirical delay data format
#'
#' @details An \code{empirical_delays} dataframe must contain (at least) two columns.
#' An 'event_date' column of type \code{Date}
#' and a 'report_delay' column of type \code{numeric}.
#' Each row represents the recording of a single delay between event and observation.
#' Typically, the 'event' here is the onset of symptoms of the disease of interest.
#' And the observation can be, for instance, case confirmation, hospital admission,
#' admission to an ICU, or death, depending on what the incidence data represents.
#' For a particular row, 'event_date' would then represent, for a single individual,
#' the date at which symptoms appeared. And 'report_delay' would represent the number
#' of time steps (as specified by \code{time_step}) until the observation was made
#' for this same individual.
#'
#' @name empirical_delay_data_format
NULL

#' Uncertainty summary
#'
#' @param uncertainty_summary_method One of these options:
#' \itemize{
#' \item{'original estimate - CI from bootstrap estimates'.
#' The confidence interval is built using bootstrapped estimates
#' and centered around the original estimates.}
#' \item{'bagged mean - CI from bootstrap estimates'.
#' The confidence interval is built using bootstrapped estimates
#' and centered around the mean of bootstrapped estimates and original estimates.}
#' }
#' @param original_values Optional. Values of reference
#' used to construct the uncertainty interval around.
#' Typically, these are estimates obtained on the original data.
#' Must be a dataframe with a timestep index column named \code{index_col}
#' and a value column named \code{value_col}.
#' The index column must not contain any \code{NA} value.
#' @param bootstrapped_values  Bootstrap
#' replicates of the original data.
#' Must be a dataframe in the long format
#' with a timestep index column named \code{index_col},
#' a bootstrap replicate index column named \code{bootstrap_id_col},
#' and a value column named \code{value_col}.
#' The index column must not contain any \code{NA} value.
#' @param central_values Values around which the confidence interval is going to be centered.
#' Must be a dataframe with a timestep index column named \code{index_col}
#' and a value column named \code{value_col}.
#' The index column must not contain any \code{NA} value.
#' @param alpha value between 0 and 1. Confidence level of the confidence interval.
#' @param combine_bootstrap_and_estimation_uncertainties boolean.
#' Combine uncertainty from Re estimation with uncertainty from observation process?
#' If \code{TRUE}, the credible intervals for Re estimates must be passed via \code{Re_HPDs}.
#' The output credible intervals
#' will be the union of bootstrapping intervals and Re estimation intervals.
#' @param Re_HPDs Optional. Credible intervals for Re estimates.
#' Use only if \code{combine_bootstrap_and_estimation_uncertainties} is \code{TRUE}.
#' @param value_col string. Name of the column containing values.
#' @param bootstrap_id_col string. Name of the column containing bootstrap samples numbering.
#' Id 0 must correspond to values associated to the original data.
#' @param index_col string. Name of the index column.
#' The index tracks which data point in bootstrapped values
#' corresponds to which data point in the original values.
#' @param output_value_col string. Name of the output column with estimated values.
#' @param prefix_up string. prefix to add to \code{output_value_col}
#' to name the column containing the upper limit of confidence interval
#' @param prefix_down string. prefix to add to \code{output_value_col}
#' to name the column containing the lower limit of confidence interval
#'
#' @name uncertainty
NULL

#' Module object utilities
#'
#' @param index_col string. Index column to keep track of which data point
#'  in bootstrapped estimates corresponds to which data point in the original estimates.
#'
#' @name module_objects
NULL
