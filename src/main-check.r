# setwd
setwd("/home/dima/pp-idirt/")

# load libs
source("src/R/exp-functions.r")
source("src/R/affinity-obj.r")
source("src/R/aff-collection-obj.r")

## load experiments

# getrate manually experiment names vector
#e.names <- c("L1_1", "L1_2", "L1_3", "1a", "1b", "1c", "2a", "2b")
# read exp template (from John)
dt.tmpl <- fread("data/tbls/experimentalDesignTemplate.txt", header = T, sep = "\t")

# get exp names
e.names <- unique(dt.tmpl$Experiment)

# load all exp in list
l.exp <- lapply(e.names, function(exp)
  mqexp(file = "data/idirt/s_total.txt", exp))
names(l.exp) <- e.names

# get affinity object for each experiment
#  (now only first protein for each protein group)
#  remove contaminants and decoy proteins
l.aff <- lapply(l.exp, function(cur.exp)
  affinity(cur.exp, minpep = 3, avalue = "H/(H+L)", use.normalized = T))

# gett affinity collection object
cur.col <- affcollection(l.aff)

# matrix with all affinity values
head(cur.col[["amtx"]])

## one experiment study

# subset one experiment from list
cur.exp <- l.exp[["1a"]]

# generate report
generate_report(cur.exp, wdir = "rep/1a")

# get affinity object
cur.aff <- affinity(cur.exp, minpep = 3, avalue = "H/(H+L)", use.normalized = T)

## get significance for each protein and plot histogram
# example without filtration
pv <- signif_test(cur.aff, apply.chauven = F,
                  add.plot = T, plot.sign = 1e-6)

# with chauvenet filtration
pv <- signif_test(cur.aff, apply.chauven = T, max.it = 10, ch.th = 2,
                  add.plot = T, plot.sign = 1e-6)
head(pv)
