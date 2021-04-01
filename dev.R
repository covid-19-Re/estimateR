library(devtools)

load_all()

devtools::document()

packagePath <- build(build_vignettes = TRUE)
install.packages(packagePath)
