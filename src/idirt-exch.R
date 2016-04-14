# set working directory
setwd("/home/dima/mega/pp-idirt")

# some libs
library(data.table)

# load setup data
source("src/idirt-setup-exch.R")
# load IDIRT functions
source("src/idirt-functions.R")

# load data from setup file
total.prj <- loadIDIRT(i.prj)

# select all avaliable experiments
work.prj <- total.prj

# find best proteins from protein groups
#v.selp <- getBestProt(work.prj)
#save(list = c("v.selp"), file = "rdat/vselp-exchange.rda")
# this function not very fast, so i save its results (for all experiments)
# to vselp.rda object and load it
load("rdat/vselp-exchange.rda")

# load vector with converstion ipi to gene names (ipi.gv)
load("rdat/ipi-gv.rda")

# get matrix with all data for selected proteins with converted names
m.plot <- getIDIRTmtx(work.prj, v.selp, ipi.gv)

# subset not na data
m.plot <- m.plot[apply(m.plot, 1, function(x) !any(is.na(x))), ]

#pdf("fig/exh_velo.pdf", width = 7, height = 9)
#par(ps = 8, mar = c(4, 4, 2, 2))
plot(x = 1:ncol(m.plot), y = m.plot[1, ], type = "l", col = "gray80",
     ylim = c(0, 1), xaxt = "n", ylab = "H / (H + L)", xlab = "",
     main = "Velos exchange")

#axis(1, at = 1:ncol(m.plot), labels = gsub(".* (.*)", "\\1", colnames(m.plot)))
for (i in 1:nrow(m.plot)) {
  lines(x = 1:ncol(m.plot), y = m.plot[i, ], col = "gray80")  
}
w.plot <- m.plot[m.plot[, 1] > 0.6, ]
#w.plot <- m.plot[sample(1:nrow(m.plot), size = 10), ]
for (i in 1:nrow(w.plot)) {
  lines(x = 1:ncol(w.plot), y = w.plot[i, ], col = "red", type = "o")  
}
text(x = 1.4, y = w.plot[, 1], labels = rownames(w.plot), col = "black",
     cex = .8)
legend("bottomleft", legend = rownames(w.plot), cex = .8)
#dev.off()

#pdf("fig/exh_velo_delta.pdf", width = 7, height = 12)
#par(ps = 12)
plot(y = m.plot[, ncol(m.plot)] - m.plot[, 2], x = m.plot[, 1], pch = 15,
     type = "n",
     xlab = colnames(m.plot)[1],
     ylab = sprintf("%s - %s delta", colnames(m.plot)[i], colnames(m.plot)[2]))
abline(h = 0, lty = 2, col = "black")

points(y = m.plot[, ncol(m.plot)] - m.plot[, 2], x = m.plot[, 1], pch = 15, cex = .9)
g.plot <- m.plot[m.plot[, 1] > 0.6, ]
text(x = g.plot[, 1], y = g.plot[, i] - g.plot[, 2] + 0.015,
     labels = rownames(g.plot),
     col = "blue", cex = .6)
#dev.off()
