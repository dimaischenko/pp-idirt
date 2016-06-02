# set working directory
setwd("/home/dima/mega/pp-idirt")

# some libs
library(data.table)

# load setup data
source("src/setup/idirt-setup.R")
# load IDIRT functions
source("src/lib/idirt-functions.R")

# load data from setup file
total.prj <- loadIDIRT(i.prj)

# load project descr
proj.desc <- fread("data/tbls/prj_info.txt", header = T, sep = "\t",
                        stringsAsFactors = F)
setkey(proj.desc, proj, exp)

# we can subset needed experiments for example:
# work.prj <- total.prj[c("IDIRT", "Tandem")]
# but now we select all avaliable experiments
work.prj <- total.prj

# find best proteins from protein groups
#v.selp <- getBestProt(work.prj)
#save(list = c("v.selp"), file = "rdat/vselp-total.rda")
# this function not very fast, so i save its results (for all experiments)
# to vselp.rda object and load it
load("rdat/vselp-total.rda")

# load vector with converstion ipi to gene names (ipi.gv)
#ipi.gv <- ipi2Gene(total.prj)
#save(list = c("ipi.gv"), file = "rdat/ipi-gv.rda")
load("rdat/ipi-gv.rda")

# get matrix with all data for selected proteins with converted names
m.plot <- getIDIRTmtx(work.prj, v.selp, ipi.gv)

## some basic plots

library(gplots)

# left proteins that has at least 0.6 value at least in one exp
m.plot <- m.plot[apply(m.plot, 1, function(x)
  any(x >= 0.6, na.rm = T)), ]
dim(m.plot)

# filter matrix to plot heatmap (we sholud exclude "NA" id dist matrix)
while(sum(is.na(dist(m.plot))) != 0) {
  tmp <- sort(apply(as.matrix(dist(m.plot)), 1, function(x) sum(is.na(x))), 
            decreasing = T)
  m.plot <- m.plot[-match(names(tmp[1]), rownames(m.plot)), ]
}
dim(m.plot)

# plot heatmap
#pdf("fig/test_hm_0.8.pdf", width = 7, height = 11)
#par(ps = 9)
heatmap.2(m.plot, Colv = F, 
          dendrogram = "row",
          trace = "none", lhei = c(2, 12),
          na.col = "white",
          breaks = seq(0, 1, length.out = 12),
          key.title = "",
          key.xlab = "H / (H + L)",
          #scale = c("none"),
          #na.rm = T,
          #symbreaks = min(m.plot, na.rm = TRUE),
          #col = c(rep("white", 11), brewer.pal(n = 11, name = "Spectral")),
          col = brewer.pal(n = 11, name = "Spectral"),
          margin = c(12, 12))
#dev.off()

## clusters

m.dist <- dist(m.plot)
hc.plot <- (hclust(dist(m.plot)))
hc.plot.cut <- cutree(hc.plot, k = 16)

#pdf("fig/all_0.8_pat.pdf", width = 8, height = 8)
#par(ps = 9)
layout(matrix(1:16, nrow = 4, byrow = T))
par(mar = c(4, 4, 2, 2))
for (k.cl in 1:16) {
plot(x = 1:ncol(m.plot), y = m.plot[1, ], type = "l",
     col = "gray60", ylim = c(0, 1), ylab = "H / (H + L)",
     xlab = "", xaxt = "n", main = sprintf("cluster: %s", k.cl))
axis(1, at = 1:ncol(m.plot), labels = rep("", ncol(m.plot)), tck = -0.035)
mtext("Experiments", side = 1, line = 1, cex = .8)
for (i in 2:nrow(m.plot)) {
  lines(x = 1:ncol(m.plot), y = m.plot[i, ], col = "gray")
}
for (k in names(hc.plot.cut[hc.plot.cut == k.cl])) {
  lines(x = 1:ncol(m.plot), y = m.plot[k, ], col = "red")
}
legend("topleft", legend = names(hc.plot.cut[hc.plot.cut == k.cl]),
       y.intersp = 1, box.col = "gray40", box.lwd = .8,
       ncol = 1, cex = .8, bg = rgb(1, 1, 1, 0.7))
}
#dev.off()
