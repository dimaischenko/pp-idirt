# Affinity collection object
# It differ from simple list of affinity objects
# because we have additional slots (e.g. affinity matrix)
# so create separate class
# Author: Dima Ischenko
# Date: 23.06.2016

#' Create affinity colelction from list of affinity objects
#' 
#' @param alist List with "affinity" objects
#' @return "affcoll" object, collection with affinit objects
#' @export
affcollection <- function(alist) {
  # create object
  acoll <- list()
  class(acoll) <- "affcoll"
  
  # add matrix
  # get all proteins
  v.prot <- unique(unlist(lapply(alist, function(x) x$protein)))
  # create empty matrix
  acoll[["amtx"]] <- matrix(NA, nrow = length(v.prot),
                                  ncol = length(alist),
                                  dimnames = list(v.prot, names(alist)))
  # fill with values
  for (i in 1:length(alist)) {
    acoll[["amtx"]][alist[[i]]$protein, names(alist)[i]] <- alist[[i]]$aff
  }
  
  return(acoll)
}