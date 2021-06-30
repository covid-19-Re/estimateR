
# Estonian case data from Terviseamet
url <- "https://opendata.digilugu.ee/opendata_covid19_tests_total.csv"
EST_data <- try(readr::read_csv(url))
if ("try-error" %in% class(EST_data)) {
  stop(stringr::str_c("Couldn't read EST case data at ", url))
}

EST_incidence_data <- EST_data %>%
  dplyr::transmute(date = lubridate::as_date(StatisticsDate),
            case_incidence = DailyCases)

usethis::use_data(EST_incidence_data, overwrite = TRUE, compress = "bzip2", version = 2)
