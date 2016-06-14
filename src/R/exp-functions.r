# New approach to experiment control
# experiments as separate objects
# Author: Dima Ischenko
# Date: 11.06.2016

# *******
# Usefull and helpfull functions
# *******

#' Simple boolean `not in` function.
#' @export
"%ni%" <- Negate("%in%")

#' Length of unique elements in vector.
#' 
#' @param x Vector with elements.
#' @return Length of unique elements.
#' @export
ulength <- function(x) length(unique(x))

#' Execute code in speciefied direcotry.
#' 
#' @param dir Directory to run code.
#' @param code Code to run.
#' @export
in_dir <- function(dir, code) {
  cur <- getwd()
  setwd(dir)
  on.exit(setwd(cur))
  
  force(code)
}

#' Chauvenet filtration of distribution.
#'
#' @param Vector with data.
#' @param Constant for normalization.
#' @return Vector without filtered values.
#' @export
ChauvenetFilter <- function(x, th = 0.5) {
  x.mean <- mean(x)
  x.sd <- sd(x)
  
  xpr <- pnorm(x, x.mean, x.sd, lower.tail = F) * length(x)
  xpl <- pnorm(x, x.mean, x.sd, lower.tail = T) * length(x)
  return(x[xpr > th & xpl > th])
}

# *******
# Main artillery
# *******

#' Load experiment from file
#' create mqexp object.
#' 
#' @param file Path to file.
#' @param exp.name Experiment name.
#' @return mqexp object.
#' @export
mqexp <- function(file, exp.name) {
  require(data.table)
  
  # load file data
  dt.a <- fread(file, header = T)
  # colnames and experiment name to upper case
  exp.Name <- toupper(exp.name)
  setnames(dt.a, toupper(colnames(dt.a)))
  # generate needed columns
  v.col <- c("MAJORITY PROTEIN IDS",
             sprintf(fmt = c("PEPTIDES %s", "RAZOR + UNIQUE PEPTIDES %s",
                             "RATIO H/L %s", "RATIO H/L NORMALIZED %s",
                             "RATIO H/L VARIABILITY [%%] %s"), exp.Name),
             "PEP", "CONTAMINANT", "REVERSE"
             )
  # report error if not all columns present in data file
  if (!all(v.col %in% colnames(dt.a))) {
    stop(sprintf("Can not load %s experiment", exp.Name))
  }
  # subset needed columns
  dt.s <- dt.a[, v.col, with = F]

  # split protein group to proteins
  # get list with splitted proteins for each protein group
  l.prot <- sapply(dt.s[["MAJORITY PROTEIN IDS"]], function(x) {
    x.u <- unlist(strsplit(x, split = ';'))
  })
  
  # generate final data
  dt.f <- dt.s[rep(1:nrow(dt.s), sapply(l.prot, length)),]
  dt.f[["PG"]] <- rep(1:nrow(dt.s), sapply(l.prot, length))
  dt.f[["PROTEIN"]] <- unlist(l.prot)
  dt.f[["MAJORITY PROTEIN IDS"]] <- NULL
  
  # rename columns
  setnames(dt.f, c("peps", "rupeps", "hl", "hln", "hlv",
                   "score", "contaminant", "reverse", "pg", "protein"))
  # manually convert needed columns to numeric
  for (col in c("peps", "rupeps", "hl", "hln", "hlv", "score")) {
    dt.f[[col]] <- as.numeric(dt.f[[col]])
  }
  
  # remove NA values
  dt.f <- dt.f[!is.na(dt.f$hl) & !is.na(dt.f$hln), ]
  
  # create and return mqexp object
  exp <- list(name = exp.name, data = dt.f)
  class(exp) <- "mqexp"
  return(exp)
}

#' Function to generate tex report for mqexp experiment.
#' 
#' @param exp mqexp object.
#' @param wdir Directory to write tex file with report.
#' @param tmpl Template for report.
#' @export
generate_report <- function(exp, wdir, tmpl = "sweave/exp-report.Rnw") {
  require(tools)
  require(utils)
  
  # get absolute path to template
  tmpl.apath <- file_path_as_absolute(tmpl)
  # create direcotry for output and generate report
  dir.create(wdir, recursive = T)
  
  # TODO(ischenko): very interesting approach
  #  but now i didn't find any solution for automatic
  #  sweave generation in new environment

  # create new driver for Sweave
  swdriver <- RweaveLatex()
  # set current environment
  swenv <- environment() 
  # redefine function to evaluation
  SwEvalWithOpt <- function (expr, options){
    if(options$eval){
      res <- try(withVisible(eval(expr, swenv)),
                 silent=TRUE)
      if(inherits(res, "try-error")) return(res)
      if(options$print | (options$term & res$visible))
        print(res$value)
    }
    return(res)
  }
  swdriver$runcode <- makeRweaveLatexCodeRunner(evalFunc = SwEvalWithOpt)
  
  # run sweave with constructed driver
  in_dir(wdir, Sweave(file = tmpl.apath, driver = swdriver))
}