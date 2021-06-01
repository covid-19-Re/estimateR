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

#' Incidence data for COVID-19 in Estonia
#'
#' A dataset containing aggregated incidence data for Estonia from February 2020 to May 2021.
#' The linelist data is published by the Estonian public health authorities.
#'
#' @format A data frame with 460 rows and 2 variables:
#' \describe{
#'   \item{report_date}{date of case reporting in YYYY-mm-dd format}
#'   \item{incidence}{number of cases reported on this date}
#' }
#' @source \url{https://opendata.digilugu.ee/opendata_covid19_tests_total.csv}
"EST_incidence_data"

#' Linelist of COVID-19 patients in Switzerland
#'
#' A dataset containing a linelist of Swiss COVID-19 patients from February to June 2020.
#'
#' Each row refers to a particular patient.
#' The 'confirmation_date' column refers to the date at which a positive COVID-19 test is reported.
#' If applicable, a date of hospital admission is included.
#' If applicable, a date of death is included.
#' If applicable and known, a date of onset of symptoms is included.
#'
#' This linelist is provided by the Federal Office of Public Health in Switzerland.
#'
#' @format A data frame with 31'950 rows and 4 columns:
#' \describe{
#'   \item{confirmation_date}{date of case confirmation in YYYY-mm-dd format}
#'   \item{symptom_onset_date}{If applicable and reported, date of onset of symptoms in YYYY-mm-dd format}
#'   \item{hospital_admission_date}{If applicable, date of admission to a hospital in YYYY-mm-dd format}
#'   \item{death_date}{If applicable, date of death of the patient in YYYY-mm-dd format}
#' }
#' @source Federal Office of Public Health
"CH_linelist"

