## code to prepare `HK_incidence_data` dataset goes here

max_date <- as.Date("2020-08-01")

HK_linelist_data_url <- "https://api.data.gov.hk/v2/filter?q=%7B%22resource%22%3A%22http%3A%2F%2Fwww.chp.gov.hk%2Ffiles%2Fmisc%2Fenhanced_sur_covid_19_eng.csv%22%2C%22section%22%3A1%2C%22format%22%3A%22csv%22%7D"
HK_linelist_data <- try(readr::read_csv(HK_linelist_data_url,
  col_types = list(
    `Case no.` = readr::col_double(),
    `Report date` = readr::col_character(),
    `Date of onset` = readr::col_character(),
    Gender = readr::col_character(),
    Age = readr::col_double(),
    `Name of hospital admitted` = readr::col_character(),
    `Hospitalised/Discharged/Deceased` = readr::col_character(),
    `HK/Non-HK resident` = readr::col_character(),
    `Case classification*` = readr::col_character(),
    `Confirmed/probable` = readr::col_character()
  )
))

if ("try-error" %in% class(HK_linelist_data)) {
  stop(stringr::str_c("Couldn't read Hong Kong linelist data at ", HK_linelist_data_url))
}

HK_linelist_data <- HK_linelist_data %>%
  dplyr::filter(.data$`Confirmed/probable` == "Confirmed") %>%
  # only keep confirmed cases
  dplyr::transmute(
    report_date = as.Date(.data$`Report date`, format = "%d/%m/%Y"),
    onset_date = as.Date(.data$`Date of onset`, format = "%d/%m/%Y")
  ) %>%
  dplyr::filter(
    !is.na(.data$report_date),
    .data$report_date < max_date
  )

# Gather the incidence based on all confirmed cases
case_incidence <- HK_linelist_data %>%
  dplyr::transmute(date = .data$report_date) %>%
  dplyr::group_by(.data$date) %>%
  dplyr::tally(name = "case_incidence")

# Gather incidence based on dates of onset of symptoms
# for cases for which this is known
onset_incidence <- HK_linelist_data %>%
  dplyr::transmute(date = .data$onset_date) %>%
  dplyr::filter(!is.na(.data$date)) %>%
  dplyr::group_by(.data$date) %>%
  dplyr::tally(name = "onset_incidence")

# Gather incidence based on dates of case confirmation
# for cases with no known date of onset of symptoms
report_incidence <- HK_linelist_data %>%
  # Only keep report_date when we do not have the onset_date
  dplyr::mutate(report_date = dplyr::if_else(is.na(.data$onset_date),
    .data$report_date,
    as.Date(NA)
  )) %>%
  dplyr::transmute(date = .data$report_date) %>%
  dplyr::filter(!is.na(.data$date)) %>%
  dplyr::group_by(.data$date) %>%
  dplyr::tally(name = "report_incidence")

HK_incidence_data <- dplyr::full_join(onset_incidence, report_incidence, by = "date") %>%
  dplyr::full_join(y = case_incidence, by = "date") %>%
  tidyr::replace_na(list(onset_incidence = 0, report_incidence = 0, case_incidence = 0)) %>%
  tidyr::complete(
    date = seq.Date(min(.data$date), # add zeroes for dates with no reported case
      max(.data$date),
      by = "days"
    ),
    fill = list(
      onset_incidence = 0,
      report_incidence = 0,
      case_incidence = 0
    )
  ) %>%
  dplyr::select(.data$date, .data$case_incidence, .data$onset_incidence, .data$report_incidence) %>%
  dplyr::arrange(.data$date)

usethis::use_data(HK_incidence_data, overwrite = TRUE, compress = "bzip2", version = 2)
