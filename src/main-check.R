# setwd
setwd("/home/dima/pp-idirt/")

# load libs
source("src/R/exp-functions.R")

for (e.name in c("1a", "1b", "1c", "2a", "2b", "2c")) {
  # load exp
  cur.exp <- idirt(file = "data/idirt/s_idirt.csv", exp.name = e.name)
  # generate report
  generate_report(odirt = cur.exp, wdir = sprintf("rep/%s", e.name))
}
