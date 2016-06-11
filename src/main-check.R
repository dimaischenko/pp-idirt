# setwd
setwd("/home/dima/pp-idirt/")

# load libs
source("src/R/exp-functions.R")

for (e.name in c("1a", "1b", "1c", "2a", "2b", "2c")) {
  # load exp
  dt.e <- LoadExp("data/idirt/s_idirt.csv", e.name)
  # generate report
  generate_report(dt.e, sprintf("rep/%s", e.name))
}