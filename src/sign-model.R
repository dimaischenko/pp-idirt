library(VGAM, quietly = T)
library(gplots, quietly = T)
library(RColorBrewer, quietly = T)

setwd("/home/dima/mega/pp-idirt/")

library(data.table)

# load setup data
source("src/idirt-setup.R")
# load IDIRT functions
source("src/idirt-functions.R")

# load data from setup file
total.prj <- loadIDIRT(i.prj)

# we can subset needed experiments for example:
work.prj <- total.prj[c("IDIRT")]

# find best proteins from protein groups
#v.selp <- getBestProt(work.prj)
#save(list = c("v.selp"), file = "rdat/vselp-total.rda")
# this function not very fast, so i save its results (for all experiments)
# to vselp.rda object and load it
load("rdat/vselp-total.rda")

# load vector with converstion ipi to gene names (ipi.gv)
load("rdat/ipi-gv.rda")

# get matrix with all data for selected proteins with converted names
m.plot <- getIDIRTmtx(work.prj, v.selp, ipi.gv)

# select two experiments
exp.pr <- c("2a", "2b")

# subset data
dt <- m.plot[, exp.pr]
# filter NA
dt <- dt[apply(dt, 1, function(x) sum(is.na(x)) == 0), ]

# number of proteins
k.Np <- nrow(dt)

# theoretical distribution
k.mu <- apply(dt, 2, mean)
k.sd <- apply(dt, 2, sd)
k.Np.th <- k.Np * 1e2 # number of points in theoretical

# create theoretical
dt.th <- do.call(cbind, lapply(1:ncol(dt), function(i)
  rnorm(k.Np.th, k.mu[i], k.sd[i])))

# calculate distance between each two proteins
x.pd <- sqrt(matrix(apply(dt, 1, function(x)
  apply(dt, 1, function(y) sum((x - y)^2))),
  nrow = k.Np))
  
# calculate p-value of the distance
x.pv <- matrix(sapply(1:nrow(dt), function(i) {
  th.ecdf <- ecdf(sqrt(apply(dt.th, 1, function(v) 
    sum((dt[i, ] - v)^2))))
  th.ecdf(x.pd[i, ])
}), nrow = k.Np)

# adjust by multiple correction
x.pv.adj <- matrix(p.adjust(x.pv, method = "BH"), nrow = k.Np)

# get significant pairs
p.sign <- 0.01 # p-value threshold
sig.pnts <- which(x.pv.adj < p.sign, arr.ind = T)
sig.pnts <- sig.pnts[sig.pnts[, 1] > sig.pnts[, 2], ]

# get all significant proteins
sgp <- unique(rownames(dt)[sig.pnts])
# plot all
plot(dt[, 1], dt[, 2], pch = 16, col = "gray", xlab = exp.pr[1],
  ylab = exp.pr[2])
# add significant
points(dt[sgp, 1], dt[sgp, 2], pch = 21, col = "black", bg = "orange")

cbind(rownames(dt)[sig.pnts[, 1]], rownames(dt)[sig.pnts[, 2]])