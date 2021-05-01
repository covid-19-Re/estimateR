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
#' @param incidence_data numeric vector. Incidence data (TODO improve definition)
#' @param simplify_output boolean. Return a numeric vector instead of module
#'  output object if output offset is zero? TODO to be described better.
#'
#' @name module_structure
NULL

#' Inner module option characteristics
#'
#' @param incidence_input module input object.
#' @name inner_module
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
#'
#' @name pipe_params
NULL

#' High-level delay parameters
#'
#' @param delay_incubation TODO add details
#' @param delay_onset_to_report TODO add details
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
#' }
#' @param uncertainty_summary_method string. One of the following options:
#' \itemize{
#' \item{'NONE' if no summary of bootstrap estimates is required}
#' \item{'original estimate - CI from bootstrap estimates'}
#' \item{'bagged mean - CI from bootstrap estimates'}
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
#' @param matrix_a,matrix_b discretized delay distribution matrix
#'
#' @name distribution
NULL

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
#' @param original_estimates Optional. TODO add details
#' @param bootstrapped_estimates TODO add details
#' @param Re_estimate_col string. Name of the column containing Re estimates
#' @param bootstrap_id_col string. Name of the column containing bootstrap samples numbering.
#' Id 0 must correspond to the estimate on the original data.
#' @param index_col string. Index column to keep track of which data point
#'  in bootstrapped estimates corresponds to which data point in the original estimates.
#' @param alpha value between 0 and 1. Confidence level.
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






