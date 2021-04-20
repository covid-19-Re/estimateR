## code to prepare `HK_incidence_data` dataset goes here

HK_linelist_data_url <- "https://api.data.gov.hk/v2/filter?q=%7B%22resource%22%3A%22http%3A%2F%2Fwww.chp.gov.hk%2Ffiles%2Fmisc%2Fenhanced_sur_covid_19_eng.csv%22%2C%22section%22%3A1%2C%22format%22%3A%22csv%22%7D"
HK_linelist_data <- try(readr::read_csv(HK_linelist_data_url,
                                        col_types = list(
                                          `Case no.` = col_double(),
                                          `Report date` = col_character(),
                                          `Date of onset` = col_character(),
                                          Gender = col_character(),
                                          Age = col_double(),
                                          `Name of hospital admitted` = col_character(),
                                          `Hospitalised/Discharged/Deceased` = col_character(),
                                          `HK/Non-HK resident` = col_character(),
                                          `Case classification*` = col_character(),
                                          `Confirmed/probable` = col_character()
                                        )))

if ("try-error" %in% class(HK_linelist_data)) {
  stop(str_c("Couldn't read Hong Kong linelist data at ", HK_linelist_data_url))
}

HK_incidence_data <- HK_linelist_data %>%
  dplyr::filter(.data$`Confirmed/probable` == "Confirmed") %>% # only keep confirmed cases
  dplyr::transmute(report_date = as.Date(.data$`Report date`, format = "%d/%m/%Y")) %>%
  # we only need the 'Report date' column (we purposefully ignore all the extra data we could get from other columns)
  dplyr::filter(!is.na(.data$report_date)) %>%
  dplyr::group_by(.data$report_date) %>%
  dplyr::tally(name = "incidence") %>% # count the incidence
  tidyr::complete(report_date = seq.Date(min(.data$report_date), # add zeroes for dates with no reported case
                                         max(.data$report_date),
                                         by = "days"),
                  fill = list(incidence = 0)) %>%
  dplyr::arrange(.data$report_date)

usethis::use_data(HK_incidence_data, overwrite = TRUE, compress = "bzip2", version = 2)
