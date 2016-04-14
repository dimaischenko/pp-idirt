# functions
ulength <- function(x) length(unique(x))
"%ni%" <- Negate("%in%")

# huge function to load data from setup list
# create matricies with data
# find IPI-gene connections
# return list in "setup list" format with additional slots
# with matricies and data.frames
loadIDIRT <- function(l.prj) {
  require(data.table)
  res.data <- l.prj
  
  # load data
  for (i in 1:length(l.prj)) {
    # load main table
    i.dt <- fread(l.prj[[i]][["path"]])
    # get number of (razor + unique) peptides for each PG
    i.dt$pep <- as.numeric(gsub("(.*?)\\;.*", "\\1",
                                i.dt[["Peptide Counts (razor+unique)"]]))
    # filter by > 2
    i.dt <- i.dt[i.dt$pep > 2, ]
    # filter by reversed proteins (decoy)
    i.dt <- i.dt[i.dt$PEP < min(i.dt$PEP[i.dt$Reverse == "+"]), ]
    
    # save main table
    res.data[[i]][["data"]] <- i.dt
    
    # work with protein names
    # split by proteins part
    v.pat <- c("^Orf", "^CON", "^REV", "^IPI")
    
    i.dt[["Majority Protein IDs"]] <- gsub(">", "", i.dt[["Majority Protein IDs"]])
    l.prot <- sapply(i.dt[["Majority Protein IDs"]], function(x) {
            x.u <- unlist(strsplit(x, split = ';'))
            v.prot <- unlist(sapply(v.pat, function(pat) {
                i.match <- grep(pat, x.u)
                gsub(sprintf("(%s.*?)\\|.*", pat), "\\1", x.u[i.match])
            }))
            })
    
    # create data for ipi to prot conversion in future
    ipi.gen <- do.call(rbind, sapply(i.dt[["Majority Protein IDs"]], function(x) {
      x.u <- unlist(strsplit(x, split = ';'))
      i.match <- grep("^IPI", x.u)
      cbind(gsub("(^IPI.*?)\\|.*", "\\1", x.u[i.match]),
        gsub("^IPI.*?\\|(.*?)[\\;|\\|].*", "\\1", x.u[i.match]))
    }, simplify = "array"))
    ipi.gen[, 2] <- gsub("(IPI.*?)\\|(.*)", "\\2", ipi.gen[, 2])
    
    # save genes names
    res.data[[i]][["ipi"]] <- ipi.gen
    
    # protein data matrix
    i.mdt <- i.dt[rep(1:nrow(i.dt), sapply(l.prot, length)),
                  l.prj[[i]][["cols"]], with = F]
    i.mdt$pg <- rep(1:nrow(i.dt), sapply(l.prot, length))
    i.mdt$prot <- unlist(l.prot)
    
    # add contaminant logical value
    i.mdt$isContam <- rep(i.dt$Contaminant == "+", sapply(l.prot, length))
    
    # reverse H/L if needed (if rev in idirt-setup for column is TRUE)
    for (kr in 1:length(l.prj[[i]][["rev"]])) {
      if (l.prj[[i]][["rev"]][kr]) {
        i.mdt[[kr]] <- 1 / i.mdt[[kr]]
      }
      i.mdt[[kr]] <- (i.mdt[[kr]]^-1 + 1)^-1
    }
    
    # save matrix
    res.data[[i]][["mdata"]] <- i.mdt
    colnames(res.data[[i]][["mdata"]])[1:length(res.data[[i]]$cols)] <- names(l.prj[[i]]$cols)
  }
  
  return(res.data)
}

# function for getting best protein for each protein group
# in current dataset. it merge all data in one data.table
# connect proteins with links if they appered in the same protein group at
# all experiments. find connected components in graph
# and subset protein with less NA values for each connected component
getBestProt <- function(l.prj) {
  require(cluster)
  require(data.table)

  # merge protein groups from all data to one data.table
  d.pg <- l.prj[[1]][["mdata"]][, c("prot", "pg"), with = F]
  if (length(l.prj) >= 2) {
    for (i in 2:length(l.prj)) {
      i.dt <- l.prj[[i]][["mdata"]][
        , c("prot", "pg"), with = F]
      setnames(i.dt, c("prot", paste("pg", i)))
      d.pg <- merge(d.pg, i.dt, by = "prot", all = T)
    }
  }

  # create matrix from data.table
  m.pg <- as.matrix(d.pg[, -1, with = F])
  rownames(m.pg) <- d.pg$prot
  
  # statistic: the number of NAs for each protein
  v.pg.na <- apply(m.pg, 1, function(x) sum(is.na(x)))

  # calculate number of links for each protein pair (when they in the same
  #  protein group)
  dst.pg <- matrix(apply(m.pg, 1, function(x) {
    apply(m.pg, 1, function(y) {
      sum(x == y, na.rm = T) / sum(!is.na(x) & !is.na(y))
    })}), nrow = nrow(m.pg),
    dimnames = list(rownames(m.pg), rownames(m.pg)))

  # prepare distance matrix (0 means protein not linked, 1 - means proteins linked)
  dst.pg[is.nan(dst.pg)] <- 1
  diag(dst.pg) <- 1
  dst.pg <- 1 - dst.pg
  dst.pg[dst.pg > 0] <- 1
  
  # create graph
  gr.pg <- graph.adjacency(dst.pg, mode = "undirected")

  # find connected components
  gr.clust <- clusters(gr.pg, mode = "strong")

  # subset best protein (with less "NA" experiments for each connected component)
  # !TODO(dima) reliable sorting to get needed names (Orf not IPI)
  v.selp <- sapply(1:gr.clust$no, function(cl) {
    na.s <- sort(v.pg.na[names(gr.clust$membership[gr.clust$membership == cl])])
    na.ms <- na.s[na.s == min(na.s)]
    na.ms[length(na.ms)]
    }
  )

  return(v.selp)
}