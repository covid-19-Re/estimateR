## code to prepare `HK_delay_data` dataset goes here

max_delay_confirm <- 30

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

HK_delay_data <- HK_linelist_data %>%
  dplyr::filter(.data$`Confirmed/probable` == "Confirmed") %>% # only keep confirmed cases
  dplyr::transmute(event_date =  as.Date(.data$`Date of onset`, format = "%d/%m/%Y"), # rename/reformat columns
                   report_date = as.Date(.data$`Report date`, format = "%d/%m/%Y")) %>%
  dplyr::filter(!is.na(.data$event_date), !is.na(.data$report_date)) %>%
  dplyr::mutate(report_delay = as.integer(.data$report_date - .data$event_date)) %>% # extract reporting delays
  dplyr::mutate(report_delay = if_else(!between(.data$report_delay, 0, max_delay_confirm), # curate negative or too large reporting delays
                         as.integer(NA),
                         .data$report_delay)) %>%
  dplyr::filter(!is.na(.data$report_delay)) %>% # remove NA values
  dplyr::select(-.data$report_date) %>% # rearrange dataset
  dplyr::arrange(.data$event_date)

usethis::use_data(HK_delay_data, overwrite = TRUE, compress = "xz", version = 2)
