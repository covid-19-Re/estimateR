#' Delay between date of onset of symptoms of COVID-19 and date of case confirmation in Hong Kong
#'
#' A dataset compiling the delay (in days) between symptom onset and case report
#' for close to 3'000 individual cases.
#' This dataset is truncated in time latest and spans from January 2020 to July 2020.
#' This data was aggregated from the publicly-available linelist data for the COVID-19 epidemic in Hong Kong.
#' The linelist data is published by the Centre for Health Protection in Hong Kong.
#'
#' @format A data frame with 2,948 rows and 2 variables:
#' \describe{
#'   \item{event_date}{date of symptom onset in YYYY-mm-dd format}
#'   \item{report_delay}{number of days between symptom onset and case report}
#' }
#' @source \url{https://www.chp.gov.hk}
"HK_delay_data"


#' Incidence data for COVID-19 in Hong Kong
#'
#' A dataset containing aggregated incidence data for Hong Kong from January 2020 to July 2021.
#' This data was  put together by aggregating the publicly-available linelist data for SARS-CoV-2 in Hong Kong.
#' The linelist data is published by the Centre for Health Protection in Hong Kong.
#'
#' @format A data frame with 196 rows and 2 variables:
#' \describe{
#'   \item{report_date}{date of case reporting in YYYY-mm-dd format}
#'   \item{incidence}{number of cases reported on this date}
#' }
#' @source \url{https://www.chp.gov.hk}
"HK_incidence_data"
