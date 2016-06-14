# setwd
setwd("/home/dima/pp-idirt/")

# load libs
source("src/R/exp-functions.r")
source("src/R/affinity-obj.r")

# load experiment
cur.exp <- mqexp(file = "data/idirt/s_24_new.csv", "14")

# generate report
generate_report(cur.exp, wdir = "rep/1a")

# get affinity object (now only first protein for each protein group)
#  remove contaminants and decoy proteins
cur.aff <- affinity(cur.exp, minpep = 3, avalue = "H/(H+L)",
                    use.normalized = T)

## get significance for each protein and plot histogram

# example without filtration
pv <- signif_test(cur.aff, apply.chauven = F,
                  add.plot = T, plot.sign = 1e-6)

# with chauvenet filtration
pv <- signif_test(cur.aff, apply.chauven = T, max.it = 10, ch.th = 2,
                  add.plot = T, plot.sign = 1e-6)
head(pv)
