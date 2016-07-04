
# methods for MethyLumiSet (methylumi)

#setClass("MethyLumiSet",
# representation("VIRTUAL", description = "character")
#)

# betaqn <- function (bn){
setMethod(
   f= "betaqn",
   signature(bn="MethyLumiSet"),
   definition=function(bn){
      history.submitted <- as.character(Sys.time())
      object <- bn
      betas(object) <- betaqn (
      bn = betas(object)
      )
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with betaqn method (wateRmelon)"
      object@history <- rbind(
         object@history,
         data.frame(
            submitted = history.submitted,
            finished = history.finished,
            command = history.command
         )
      )
      object
   }
)

# naten <- function (mn, un, fudge=100 ) {  # beta1
setMethod(
   f= "naten",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100){
      history.submitted <- as.character(Sys.time())   
         object <- mn
         norm <- naten (
         mn = methylated(object), 
         un = unmethylated(object),
         fudge,
         ret2=TRUE
      ) 
      betas(object) <- norm$beta 
      methylated(object)   <- norm$methylated
      unmethylated(object) <- norm$unmethylated
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with naten method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
      object
   }
)

#nanet <- function(mn, un, fudge=100){ dyebuy1
setMethod(
   f= "nanet",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100){
      history.submitted <- as.character(Sys.time())   
         object <- mn
         norm <- nanet (
         mn = methylated(object), 
         un = unmethylated(object),
         fudge, ret2=TRUE
       ) 
      betas(object)        <- norm$beta
      methylated(object)   <- norm$methylated
      unmethylated(object) <- norm$unmethylated
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with nanet method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
      object
   }
)

fot <- function(x){
   ds <- grep( 'DESIGN', colnames(x@featureData@data), ignore.case = TRUE  )
#   stopifnot ( length(ds) == 1 )
#   ds
# Seeing as I accidently introduced this bug, until I correct the manifest
# this is my solution.
    ds[1]
}



#nanes <- function(mns, uns, onetwo, fudge=100){ dyebuy2.R
setMethod(
   f= "nanes",
   signature(mns="MethyLumiSet"),
   definition=function(mns, fudge=100){
      history.submitted <- as.character(Sys.time())   
         object <- mns
         ds <- fot(mns)
         norm <- nanes (
         mns    = methylated(object), 
         uns    = unmethylated(object),
         onetwo = object@featureData@data[,ds],
         fudge, ret2=TRUE
       ) 
      betas(object)        <- norm$beta
      methylated(object)   <- norm$methylated
      unmethylated(object) <- norm$unmethylated
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with nanes method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
      object
   }
)

#danes <- function(mn, un, onetwo, fudge=100, ...){ # dyebuy3.R
setMethod(
   f= "danes",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100, ...){
      history.submitted <- as.character(Sys.time())   
         object <- mn
         mn <- methylated(object)
         ds <- fot(object)
         norm <- danes (
         mn     , 
         un     = unmethylated(object),
         onetwo = object@featureData@data[,ds],  
         fudge, ret2=TRUE, ...
       ) 
      betas(object)        <- norm$beta
      methylated(object)   <- norm$methylated
      unmethylated(object) <- norm$unmethylated
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with danes method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
      object
   }
)

#danet <- function(mn, un, onetwo, fudge=100, ...){ # dyebuy4.R
setMethod(
   f= "danet",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100, ...){
      history.submitted <- as.character(Sys.time())   
      object <- mn
      ds <- fot(mn)
      mn <- methylated(object) 
      norm <- danet (
         mn , 
         un     = unmethylated(object),
         onetwo = object@featureData@data[,ds],
         fudge, ret2=TRUE, ...
      )
      betas(object)        <- norm$beta
      methylated(object)   <- norm$methylated
      unmethylated(object) <- norm$unmethylated
      history.finished <- as.character(Sys.time())
      history.command  <- "Normalized with danet method (wateRmelon)"
      object@history   <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished  = history.finished, 
            command   = history.command
         )
      )
      object
   }
)

#daten1 <- function(mn, un, onetwo, fudge=100, ...){ # bgeqqn
setMethod(
   f= "daten1",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100, ...){
   history.submitted <- as.character(Sys.time())   
         object <- mn
         ds     <- fot(mn)
         mn     <- methylated(object) 
         norm   <- daten1 (
            mn , 
            un     = unmethylated(object),
            onetwo = object@featureData@data[,ds],
            fudge, ret2=TRUE, ...
       ) 
      betas(object)        <- norm$beta
      methylated(object)   <- norm$methylated
      unmethylated(object) <- norm$unmethylated
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with daten1 method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished  = history.finished, 
            command   = history.command
         )
      )
      object
   }
)

#daten2 <- function(mn, un, onetwo, fudge=100, ...){  # bgeqq2.R
setMethod(
   f= "daten2",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100, ...){
      history.submitted <- as.character(Sys.time())   
      object <- mn
      ds <- fot(mn)
      mn <- methylated(object) 
      norm <- daten2 (
         mn , 
         un = unmethylated(object),
         onetwo = object@featureData@data[,ds],
         fudge, ret2=TRUE, ...
     )
     betas(object)        <- norm$beta
     methylated(object)   <- norm$methylated
     unmethylated(object) <- norm$unmethylated
     history.finished <- as.character(Sys.time())
      history.command <- "Normalized with daten2 method (wateRmelon)"
      object@history  <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished  = history.finished, 
            command   = history.command
         )
      )
      object
   }
)

#ot <- function ( mns, uns, onetwo, fudge=100) {
setMethod(
   f= "nasen",
   signature(mns="MethyLumiSet"),
   definition=function(mns, fudge=100){
         history.submitted <- as.character(Sys.time())   
         object <- mns
         ds <- fot(mns)
         norm <- nasen (
            mns = methylated(object), 
            uns = unmethylated(object),
            onetwo = object@featureData@data[,ds],
            fudge, ret2=TRUE
         )
         betas(object)        <- norm$beta
         methylated(object)   <- norm$methylated
         unmethylated(object) <- norm$unmethylated
         history.finished <- as.character(Sys.time())
         history.command <- "Normalized with nasen method (wateRmelon)"
         object@history <- rbind(
           object@history, 
           data.frame(
              submitted = history.submitted, 
              finished = history.finished, 
              command = history.command
           )
        )
        object
   }
)

#dasen <- function(mns, uns, onetwo, fudge=100, ...){
setMethod(
   f= "dasen",
   signature(mns="MethyLumiSet"),
   definition=function(mns, fudge=100, roco=NULL){
         history.submitted <- as.character(Sys.time())   
         object <- mns
         ds <- fot(mns)
         norm <- dasen (
            mns = methylated(object), 
            uns = unmethylated(object),
            onetwo=mns@featureData@data[,ds],
            fudge, roco, ret2=TRUE 
         ) 
         betas(object)        <- norm$beta
         methylated(object)   <- norm$methylated
         unmethylated(object) <- norm$unmethylated
         history.finished <- as.character(Sys.time())
         history.command <- "Normalized with dasen method (wateRmelon)"
         object@history <- rbind(
            object@history, 
            data.frame(
               submitted = history.submitted, 
               finished = history.finished, 
               command = history.command
            )
         )
         object
   }
)

#danen <- function ( mns, uns, onetwo, fudge=100, ...){
setMethod(
   f= "danen",
   signature(mns="MethyLumiSet"),
   definition=function(mns, fudge=100, ...){
   history.submitted <- as.character(Sys.time())   
         object <- mns
         ds <- fot(mns)
         norm <- danen (
         mns = methylated(object), 
         uns = unmethylated(object),
         onetwo=mns@featureData@data[,ds],
         fudge, ret2=TRUE, ... 
       ) 
         betas(object)        <- norm$beta
         methylated(object)   <- norm$methylated
         unmethylated(object) <- norm$unmethylated
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with danen method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
      object
   }
)

#tost <- function( mn, un, da, pn ) {  #pasteque
setMethod(
   f= "tost",
   signature(mn="MethyLumiSet"),
   definition=function(mn){
      history.submitted <- as.character(Sys.time())   
         object <- mn
         mn <- methylated(object) 
         betas(object) <- tost (
            mn, 
            un = unmethylated(object),
 #           da = object@featureData@data,
	    da = dfort(featureNames(object)),
            pn = pvals(object)
         ) 
#      unmethylated(object) <- NULL
#      methylated(object) <- NULL
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with tost method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
      object
   }
)

#peak.correction <- function (data, anno) {
setMethod(
   f= "fuks",
   signature(data="MethyLumiSet"),
   definition=function(data){
   history.submitted <- as.character(Sys.time())   
         object <- data
         betas(object) <- fuks (
         data = betas(object), 
         anno = object@featureData@data
       ) 
  #    unmethylated(object) <- NULL
  #    methylated(object) <- NULL
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with fuks method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
      object
   }
)

#swan <- function (mn, un, qc ) {  swan2.R
setMethod(
   f= "swan",
   signature(mn="MethyLumiSet"),
   definition=function(mn, da=NULL){
      history.submitted <- as.character(Sys.time())   
      object <- mn
      mn <- methylated(object)
      norm <- swan (
         mn,
         un = unmethylated(object), 
         qc = intensitiesByChannel(QCdata(object)),
         da,
         return.MethylSet=TRUE
       )
      methylumi:::betas(object) <- getBeta(norm)
      methylumi:::methylated(object) <- getMeth(norm)
      methylumi:::unmethylated(object) <- getUnmeth(norm)
      fData(object) <- fData(object)[rownames(betas(object)),] 
      history.finished <- as.character(Sys.time())
      history.command <- "Normalized with swan method (wateRmelon)"
      object@history <- rbind(
         object@history, 
         data.frame(
            submitted = history.submitted, 
            finished = history.finished, 
            command = history.command
         )
      )
      object
   }
)

#colnames <- function (x, do.NULL = TRUE, prefix = "col") 
setMethod(
   f= "colnames",
   signature(x="MethyLumiSet"),
   definition=function(x, do.NULL=TRUE, prefix=NULL){
      colnames(methylated(x))
   }
)



#genki <- function(bn, g=getsnp(rownames(bn)), se=TRUE ){
setMethod(
   f= "genki",
   signature(bn="MethyLumiSet"),
   definition=function(bn, se=TRUE){
      object <- bn
      bn     <- betas(object)
      g      <- getsnp(rownames(bn))
      genki( bn, g, se ) 
     
   }
)

#dmrse <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse",
   signature(betas="MethyLumiSet"),
   definition=function(betas, idmr=iDMR()){
      object    <- betas
      betas     <- betas(object)
      dmrse( betas, idmr ) 
     
   }
)

#dmrse_row <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse_row",
   signature(betas="MethyLumiSet"),
   definition=function(betas, idmr=iDMR()){
      object    <- betas
      betas     <- betas(object)
      dmrse_row( betas, idmr ) 
     
   }
)

#dmrse_col <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse_col",
   signature(betas="MethyLumiSet"),
   definition=function(betas, idmr=iDMR()){
      object    <- betas
      betas     <- betas(object)
      dmrse_col( betas, idmr ) 
     
   }
)

#seabi <- function (bn, stop, sex, X){
setMethod(
   f= "seabi",
   signature(bn="MethyLumiSet"),
   definition=function( bn, stop=1, sex, X ){
      object    <- bn
      bn     <- betas(object)
      seabi( bn, stop, sex, X ) 
     
   }
)
#pfilter<-function(mn, un, bn, da, onetwo, pn, bc, perCount, pnthresh, perc, pthresh){   
# filter function by Ruth Pidsley

setMethod(
   f= "pfilter",
   signature(mn="MethyLumiSet"),
   definition=function( mn,
 perCount = NULL, pnthresh = NULL, perc = NULL, pthresh = NULL ){
      
      object <- mn
      bn     <- betas(object)
#      bc     <- betas(object)
      if(exists("NBeads", assayData(object))){
        bc       <- assayData(object)$NBeads
        bc[bc<3] <- NA
      } else { 
        bc       <- betas(object)
      }
      mn     <- methylated(object)
      un     <- unmethylated(object)
      pn     <- pvals(object)
      da     <- object@featureData@data    
      l      <- pfilter (
         mn=mn, un=un, bn=bn, da=da, 
         pn=pn, bc=bc, perCount, pnthresh,perc,
         pthresh, logical.return=TRUE
      ) 
      object <- object[, l$samples]
      object[l$probes,]
   }
)


#BMIQ <- function(beta.v,design.v,nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1){
setMethod(
   f= "BMIQ",
   signature(beta.v="MethyLumiSet"),
   definition=function(
      beta.v,
      nL=3, doH=TRUE, nfit=5000,
      th1.v=c(0.2,0.75), th2.v=NULL, 
      niter=5, tol=0.001, plots=FALSE, 
      pri=FALSE 
   ){
         history.submitted <- as.character(Sys.time())
         object <- beta.v
         ds <- fot(object)
         d <- as.numeric(factor(object@featureData@data[,ds]))
         ibetas <- betas(object)
         betas <- sapply ( 
            colnames(ibetas), 
            function(name){
               ou <- try(
                  BMIQ(
                     ibetas[,name],
                     design.v=d, nL, doH, 
                     nfit, th1.v, th2.v, 
                     niter, tol, plots, 
                     sampleID=name, 
                     pri=FALSE 
                  )
               )
               if(inherits(ou, 'try-error')){
                  ou <- rep(NA,dim(ibetas)[1])
                  warn(paste(name, "failed!"))
               }else{
               ou <- ou$nbeta
               }
             }
         )
         #browser()
         colnames(betas)  <- colnames(betas(object))
         betas(object)    <- betas
         history.finished <- as.character(Sys.time())
         history.command <- "Betas processed BMIQ method "
         object@history <- rbind(
            object@history,
            data.frame(
               submitted = history.submitted,
               finished = history.finished,
               command = history.command
            )
         )
         object
   }
)



# bscon <- function(x, ...) { UseMethod (bscon, x ) }
# see also bscon_methy and bscon_minfi

setMethod(
   f= "bscon",
   signature(x="MethyLumiSet"),
   definition=function( x ){
      bscon_methy(x)
   }
)

# outlyx <- function(x, full, iqr, iqrP, pc, mv, mvP)
setMethod(
   f= "outlyx",
   signature(x="MethyLumiSet"),
   definition=function(x, iqr, iqrP, pc, mv, mvP, plot){
   x <- betas(x)
   outlyx(x, iqr, iqrP, pc, mv, mvP, plot)
   }
)

setMethod(
   f= "pwod",
   signature(object="MethyLumiSet"),
   definition=function(object, mul){
     history.submitted <- as.character(Sys.time())
     object <- betas(object)
     newbetas <- pwod(object, mul)
     betas(object) <- newbetas
     history.finished <- as.character(Sys.time())
     history.command <- "Betas processed with pwod"
     object@history <- rbind(
       object@history,
       data.frame(
         submitted = history.submitted,
         finished = history.finished,
         command = history.command
         )
       )
   object
   }
)

setMethod(
  f= "agep",
  signature(betas="MethyLumiSet"),
  definition=function(betas, coeff = NULL, verbose = FALSE){
    object <- betas(betas)
    agep(betas=object, coeff, verbose)
  }
)
   
