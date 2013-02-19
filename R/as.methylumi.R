# generic takes naked data. gets annotation from IlluminaHumanMethylation450k.db
# takes retains whatever data available (does not recalculate betas)
# methylumiset method would be useful to complete annotation

as.methylumi <- function( 
   mn=NULL, 
   un=NULL, 
   bn=NULL, 
   pv=NULL, 
   qc=NULL, 
   da=NULL,
   fd=c('CHR','DESIGN')  
) {
   dt <-  list ( mn, un, bn, pv, da )
   data <- ! sapply ( dt, is.null )
   if (! any (data) ) {
      stop('no data!')}
   rn <- rownames(dt[[ which(data)[1] ]])
   if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
      stop('can\'t load methylumi package')}
   history.submitted <- as.character(Sys.time())   
   x <- new("MethyLumiSet")
   if(!is.null(mn))methylated(x)      <- mn
   if(!is.null(pv))pvals(x)           <- pv
   if(!is.null(bn))betas(x)           <- bn
   if(!is.null(un))unmethylated(x)    <- un
   if(!is.null(qc))QCdata(x)          <- qc
   if (is.null(da))da <- data.frame(pop(fd, rn))
   fData(x) <- da
      history.finished <- as.character(Sys.time())
      history.command <- "created with as.methylumi (wateRmelon)"
      x@history <- rbind(
         x@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
   x  
}

pop <- function (fd, rn){
   o <- data.frame(row.names=rn)
   for (col in fd){
      thing  <- paste("IlluminaHumanMethylation450k", col, sep='')
      stuff <- get(thing)
      log  <- rn %in% keys(stuff)
      data <- toTable(stuff[rn[log]])
      #for (item in colnames(data)[-1]){
      #   o[log,item] <- data[,item]
      #}
      o[data[,1],colnames(data)[-1]] <- data[,-1]
   }
   o
}



setGeneric ( name= "as.methylumi"    )

setMethod(
   f= "as.methylumi",
   signature(mn="MethyLumiSet"),
   definition= function( 
      mn,
      un=NULL, 
      bn=NULL, 
      pv=NULL, 
      qc=NULL, 
      da=NULL,
      fd=c('CHR','DESIGN')  
   ) {
  object <- mn
  as.methylumi(
      mn=methylated(object),
      un=unmethylated(object), 
      bn=betas(object), 
      pv=pvals(object), 
      qc=QCdata(object), 
      da=NULL,
      fd=c('CHR','DESIGN')  
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
      da=NULL,
      fd=c('CHR','DESIGN')  
   ) {
  object <- mn
  mn <- getMeth(object)
  un <- getUnmeth(object) 
  bn <- getBeta(object) 
  #browser()
  as.methylumi(
      mn,
      un, 
      bn, 
      pv=NULL,   # need RGSet for this 
      qc=NULL,   # probly need RGSet for this
      da=NULL,
      fd=c('CHR','DESIGN')  
   )
}   
)


getColumns <- function(){
   gsub(
      'IlluminaHumanMethylation450k', 
      '', 
      ls("package:IlluminaHumanMethylation450k.db")
   )
}



#   if (!is.null(pv)) {
#      if (!all.equal(rownames(pv), rownames(methylated)) {
#        stop (pv not in the same order as mn)
#      }
#   } 
