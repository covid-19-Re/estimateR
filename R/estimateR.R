#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom rlang :=
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
#' @param uncertainty_summary_method string. TODO continue here
#' \itemize{
#' \item{'original estimate - CI from bootstrap estimates'}
#' \item{'bagged mean - CI from bootstrap estimates'}
#' } 'NONE' if no summary of bootstrap estimates is made
#' or one of the possible strings in \code{\link{summarise_uncertainty}}
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
#'
#' @name distribution
NULL






