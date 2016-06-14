# setwd
setwd("/home/dima/pp-idirt/")

# load libs
source("src/R/exp-functions.R")

# load experiment
cur.exp <- idirt(file = "data/idirt/s_24_new.csv", "12")

# generate report
generate_report(cur.exp, wdir = "rep/1a")