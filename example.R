
rm(list = ls())

packages <- c("dplyr", "vars", "tibble", "broom", "matrixStats")

for ( i in 1:length(packages) ) {
  eval(bquote(if(!require(.(packages[i]))){install.packages(.(packages[i]), dependencies = TRUE)}))
  library(packages[i], character.only = TRUE)
}



source('lag_maker.R')
source('COVIDml_rho2.R')
source('COVID_var.R')


data <- data.frame(Date = seq.Date(from = as.Date("2014-01-01"), by = "months", length.out = 100),
                   x1 = rnorm(100),
                   x2 = rnorm(100))


var <- covid_var(data, 2, "const", 50:51)
