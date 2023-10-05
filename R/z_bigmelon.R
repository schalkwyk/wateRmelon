
#{{{   epicv2clean S3 default method

#' Strip and subset EPICv2 data to work with legacy data and methods
#' 
#' @description
#' Returns an object with rownames stripped of the EPICv2 suffixes, duplicate probes are omitted.
#'
#' @details
#' EPICv2 manifests contain a few thousand probes with up to 10 replicate syntheses. 
#' To accomodate this a modified naming scheme is used, so none of the probe names match
#' those on the EPIC and previous arrays (even though most of the probes are the same sequence
#' and presumably simiar performance). 
#' 
#' This simple function relies on the rowname and subsetting methods and will work for matrix, 
#' dataframe, MethyLumiSet, or MethylSet objects, and there is a method for gds (bigmelon) objects.

epicv2clean.default <- function(x){ 
  # allow use of EPIC/450K/27K coefficients with EPICv2 data by stripping suffix
  # and discarding duplicated probenames. Names in ageCoefs don't include any of the dups.
  # x can be anything subsettable with rownames
    
    rn <- rownames(x)
    rn <- gsub('_.*$','',rn)
    go <- !duplicated(rn)
    x  <- x[go,]
    rownames(x) <- rn[go]
    x
}

# }}}
#{{{ epicv2clean S3 method for gds

epicv2clean.gds.class <- function (x) 
{
    rn <- rownames(x)
    rn <- gsub("_.*$", "", rn)
    go <- !duplicated(rn)
    
   chainsaw(gfile=x,i=go,j='')

   # x <- x[go, ]
   # rownames(x) <- rn[go]
   # x
}
# }}}
