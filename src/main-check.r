# setwd
setwd("/home/dima/pp-idirt/")

# load libs
source("src/R/exp-functions.r")
source("src/R/affinity-obj.r")

# load experiment
cur.exp <- mqexp(file = "data/idirt/s_24_new.csv", "12")

# generate report
generate_report(cur.exp, wdir = "rep/1a")

# get affinity object
cur.aff <- affinity(cur.exp, minpep = 3, avalue = "H/(H+L)",
                    use.normalized = T)