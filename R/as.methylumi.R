
as.methylumi <- function( 
   mn=NULL, 
   un=NULL, 
   bn=NULL, 
   pv=NULL, 
   qc=NULL, 
   da=NULL, ...
) {
   dt <-  list ( mn, un, bn, pv, qc, da )
   data <- ! sapply ( dt, is.null )
   if (! any (data) ) {
      stop('no data!')
   }
   history.submitted <- as.character(Sys.time())
   y <- assayDataNew( 'environment',
     betas = bn
   )
   if(!is.null(mn))y$methylated      <- mn
   if(!is.null(un))y$unmethylated    <- un
   if(!is.null(pv))y$pvals           <- pv
   if(!is.null(qc))y$QCdata          <- qc
   x <- new("MethyLumiSet", assayData=y, fData=da)
   fData(x) <- da
   history.finished <- as.character(Sys.time())
   history.command <- "created with as.methylumi (wateRmelon)"
   x@history <- rbind(
     x@history, 
     data.frame(
       submitted = history.submitted, 
       finished  = history.finished, 
       command   = history.command
     )
   )
#browser()
   x  
}


setMethod(
   f= "as.methylumi",
   signature(mn="MethyLumiSet"),
   definition= function( 
      mn,
      un=NULL, 
      bn=NULL, 
      pv=NULL, 
      qc=NULL, 
      da=NULL
   ) {
  object <- mn
  as.methylumi(
      mn=methylated(object),
      un=unmethylated(object), 
      bn=betas(object), 
      pv=pvals(object), 
      qc=QCdata(object), 
      da=NULL
   )   
   }
)


setMethod(
   f= "as.methylumi",
   signature(mn="MethylSet"),
   definition= function( 
      mn,
      un=NULL, 
      bn=NULL, 
      pv=NULL, 
      qc=NULL, 
      da=NULL
   ) {
  object <- mn
  mn  <- getMeth(object)
  un  <- getUnmeth(object) 
  bn  <- getBeta(object) 
  ann <-data.frame(getAnnotation(object), stringsAsFactors=F)
  ann$DESIGN <- ann$Type
  as.methylumi(
      mn,
      un, 
      bn, 
      pv=NULL,
      qc=NULL,
      da=ann
   )
#   browser()
  }   
)



# get the current Illumina annotation file 


aoget <- function(url= paste(
   'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/',
   'HumanMethylation450_15017482_v1-2.csv', sep='')
   ) {
   op <- options(stringsAsFactors=FALSE)
   ao <- read.csv(url, skip =7 )
   options(op)
   
}
