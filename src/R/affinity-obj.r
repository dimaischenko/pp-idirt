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
#' @param normalization Type of normalization (now not used)
#' @export
affinity.mqexp <- function(exp, minpep = 3, avalue = "H/(H+L)",
                           use.normalized = T,
                           normalization = "") {
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
  dt$val <- ifelse(rep(use.normalized, nrow(dt)),
                   dt$hln, dt$hl)
  # calculate affinity value
  dt$aff <- (dt$val^-1 + 1)^-1
  return(dt[, c("protein", "pg", "aff"), with = F])
}

#' @export
affinity.default <- function(exp, ...) {
  stop(sprintf("This class (%s) of experiment doesn't supported.",
       class(exp)))
}

