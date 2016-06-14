# Affinity class introduced to work with 
# cleaned MaxQuant results, wihtout contaminants
# normalized, converted H/L to needed value e.t.c.
# Author: Dima Ischenko
# Date: 14.06.2016

#' Convert ms experiment object to affinity object
#' 
#' @param exp MS experiment object
#' @note Now supported only "mqexp" objects
#' 
#' @return Affinity object
#' @export
affinity <- function(exp, ...) UseMethod("affinity")

#'
#' @param minpep Minimal number of peptides to consider protein
#' @param avalue Which value use to calculate affinity:
#'   H/L, L/H, H/(H+L) or L/(H+L) (!not used now)
#' @param use.normalized Use or not Normalized ratios of H/L to
#'   calculate affinity
#' @param protnames Specified names of selected proteins from
#'   mqexp object, if protnames is empty character ("") subset
#'   first protein for each protein group
#' @export
affinity.mqexp <- function(exp, minpep, avalue,
                           use.normalized,
                           protnames = c(""), ...) {
  # check avalue
  posv <- c("H/L", "L/H", "H/(H+L)", "L/(H+L)")
  if (! avalue %in% posv) {
    stop("Wrong avalue type, possible: H/L, L/H, H/(H+L) or L/(H+L)")
  }
  # get data
  dt <- exp[["data"]]
  # calc minimal decoy p-score
  mn.pep <- min(dt[["score"]][dt[["reverse"]] == "+"])
  # filter experiment data
  dt <- dt[dt$score < mn.pep & dt$peps >= minpep & 
             dt$contaminant == "", ]
  # prepare affinity
  dt$val <- ifelse(rep(use.normalized, nrow(dt)),
                   dt$hln, dt$hl)
  if (avalue %in% c("L/H", "L/(H+L)")) dt$val <- 1 / dt$val
  if (avalue %in% c("H/(H+L)", "L/(H+L)")) dt$val <- (dt$val^-1 + 1)^-1
  dt$aff <- dt$val
  
  # subset proteins
  if (protnames == "") {
    dt <- dt[, .SD[1], by = list(pg)]
  } else {
    dt <- dt[dt$protein %in% protnames, ]
  }
  aff <- dt[, c("protein", "aff"), with = F]
  class(aff) <- c("affinity", class(aff))
  
  return(aff)
}

#' @export
affinity.default <- function(exp, ...) {
  stop(sprintf("This class (%s) of experiment doesn't supported.",
       class(exp)))
}

#' Check for significance for each protein in affinity object
#'
#' @param dt Affinity object
#' 
#' @export
signif_test <- function(dt, ...) UseMethod("signif_test")

#' 
#' @param apply.chauven Logical, apply or not Chauvenet filtration for
#'   affinity value vector
#' @param maxit Maximal number of iteration in Chauvenet filtration procedure
#' @param ch.th Chauvenet filtration coefficient
#' @param add.plot Plot histogram with significance
#' @param plot.sign Significance to mark values at histogram
#' 
#' @return vector with significance for each protein
#' @export
signif_test.affinity <- function(dt, apply.chauven, max.it = 10, ch.th = 2,
                                 add.plot, plot.sign = 1e-2, ...) {
  tval <- dt$aff
  
  # Chauvenet filtration
  if (apply.chauven) {
    n.it <- 0
    while(n.it < max.it) {
      tval.n <- ChauvenetFilter(tval, th = ch.th)
      if (length(tval.n) == length(tval)) break
      tval <- tval.n
      n.it <- n.it + 1
    }
  }
  dt.mean <- mean(tval)
  dt.sd <- sd(tval)
  pv <- pnorm(q = dt$aff, mean = dt.mean, sd = dt.sd, lower.tail = F)
  names(pv) <- dt$protein
  
  # add histogram
  pv.th <- plot.sign
  hl.th <- qnorm(1 - pv.th, mean = dt.mean, sd = dt.sd)
  if (add.plot) {
    v.brks <- seq(min(dt$aff), max(dt$aff), length.out = 50)
    v.clr <- rep("white", 50)
    v.fdb <- min(which(v.brks >= hl.th))
    if (v.fdb == Inf) { v.fdb <- length(v.brks)}
    v.clr[(v.fdb - 1):length(v.clr)] <- "dodgerblue"
    hist(dt$aff, prob = T,
         main = "Affinity histogram", xlab = "Affinity",
         breaks = v.brks, col = v.clr)
    lines(density(rnorm(n = 10*nrow(dt), mean = dt.mean, sd = dt.sd)),
          col = "red", lwd = 2)
    abline(v = hl.th, lty = 2, col = "gray30")
  }
  
  return(pv)
}

#' @export
signif_test.default <- function(dt, ...) {
  stop(sprintf("This class (%s) of experiment doesn't supported.",
               class(dt)))
}