#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom rlang :=
NULL

#' Utility functions for input validity.
#'
#' @param string_user_input A string containing the value that the user passed for the tested string type parameter.
#' @param parameter_name A string containing the name the tested parameter had in the initial function in which it was passed.
#' @param incidence_data_length A number representing the length of the given incidence data.
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
#' @param partially_delayed_incidence TODO add details
#' @param fully_delayed_incidence TODO add details
#' @param simplify_output boolean. Return a numeric vector instead of module
#'  output object if output offset is zero?
#'
#'  @return A list with two elements:
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

#' Inner module option characteristics
#'
#' @param incidence_input,input,output Module input object.
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
#' @param mean_serial_interval Numeric positive value. \code{mean_si} for \code{\link[EpiEstim]{estimate_R()}}
#' @param std_serial_interval Numeric positive value. \code{std_si} for \code{\link[EpiEstim]{estimate_R()}}
#' @param mean_Re_prior Numeric positive value. \code{mean prior} for \code{\link[EpiEstim]{estimate_R()}}
#'
#' @return If \code{output_HPD = FALSE},
#' value is a module object (a list of the same kind as \code{incidence_input}).
#' The \code{values} element of the list then contains the Re estimates.
#' If \code{output_HPD = TRUE}, a list of three module objects is returned.
#' \itemize{
#' \item \code{Re_estimate} contains the Re estimates.
#' \item \code{Re_highHPD} and \code{Re_lowHPD} contain
#' the higher and lower boundaries of the HPD interval,
#' as computed by \code{\link[EpiEstim]{estimate_R()}}
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
#' @param output_infection_incidence_only boolean. TODO add details
#'
#' @name pipe_params
NULL

#' Bootstrapping parameters
#'
#' @param N_bootstrap_replicates integer. Number of bootstrap samples.
#'
#' @name bootstrap_params
NULL

#' High-level delay parameters
#'
#' @param delay_incubation TODO add details
#' @param delay_onset_to_report TODO add details
#' @param delay_distribution_final_report TODO add details
#' @param is_partially_reported_data boolean TODO add details
#' @param delay_until_partial TODO add details
#' @param delay_from_partial_to_full TODO add details
#' @param partial_observation_requires_full_observation boolean TODO add details
#'
#' @name delay_high
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
#' @param value_col string. Name of the column containing values.
#' @param bootstrap_id_col string. Name of the column containing bootstrap samples numbering.
#' Id 0 must correspond to values associated to the original data.
#' @param index_col string. Name of the index column.
#' The index tracks which data point in bootstrapped values
#' corresponds to which data point in the original values.
#' @param alpha value between 0 and 1. Confidence level of the confidence interval.
#' @param combine_bootstrap_and_estimation_uncertainties boolean.
#' Combine uncertainty from Re estimation with uncertainty from observation process?
#' If \code{TRUE}, the credible intervals for Re estimates must be passed via \code{Re_HPDs}.
#' The output credible intervals
#' will be the union of bootstrapping intervals and Re estimation intervals.
#' @param Re_HPDs Optional. Credible intervals for Re estimates.
#' Use only if \code{combine_bootstrap_and_estimation_uncertainties} is \code{TRUE}.
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
