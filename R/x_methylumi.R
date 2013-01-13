
# methods for MethyLumiSet (methylumi)

#setClass("MethyLumiSet",
# representation("VIRTUAL", description = "character")
#)


# betaqn <- function (bn){ 
setMethod(
   f= "betaqn",
   signature(bn="MethyLumiSet"),
   definition=function(bn){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
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

# naten <- function (mn, un, fudge=100 ) {
setMethod(
   f= "naten",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      history.submitted <- as.character(Sys.time())   
         object <- mn
         betas(object) <- naten (
         mn = methylated(object), 
         un = unmethylated(object),
         fudge
       ) 
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

#nanet <- function(mn, un, fudge=100){
setMethod(
   f= "nanet",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      history.submitted <- as.character(Sys.time())   
         object <- mn
         betas(object) <- nanet (
         mn = methylated(object), 
         un = unmethylated(object),
         fudge
       ) 
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
   ds <- grep( 'DESIGN', colnames(x@featureData@data) )
   stopifnot ( length(ds) == 1 )
   ds
}



#nanes <- function(mns, uns, onetwo, fudge=100){
setMethod(
   f= "nanes",
   signature(mns="MethyLumiSet"),
   definition=function(mns, fudge=100){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      history.submitted <- as.character(Sys.time())   
         object <- mns
         ds <- fot(mns)
         betas(object) <- nanes (
         mns    = methylated(object), 
         uns    = unmethylated(object),
         onetwo = object@featureData@data[,ds],
         fudge
       ) 
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

#danes <- function(mn, un, onetwo, fudge=100, ...){ 
setMethod(
   f= "danes",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100, ...){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      history.submitted <- as.character(Sys.time())   
         object <- mn
         mn <- methylated(object)
         ds <- fot(object)
         betas(object) <- danes (
         mn     , 
         un     = unmethylated(object),
         onetwo = object@featureData@data[,ds],  
         fudge, ...
       ) 
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

#danet <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "danet",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100, ...){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      history.submitted <- as.character(Sys.time())   
      object <- mn
      ds <- fot(mn)
      mn <- methylated(object) 
      betas(object) <- danet (
         mn , 
         un     = unmethylated(object),
         onetwo = object@featureData@data[,ds],
         fudge, ...
       ) 
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

#daten1 <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "daten1",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100, ...){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
   history.submitted <- as.character(Sys.time())   
         object <- mn
         ds     <- fot(mn)
         mn     <- methylated(object) 
         betas(object) <- daten1 (
            mn , 
            un     = unmethylated(object),
            onetwo = object@featureData@data[,ds],
            fudge, ...
       ) 
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

#daten2 <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "daten2",
   signature(mn="MethyLumiSet"),
   definition=function(mn, fudge=100, ...){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      history.submitted <- as.character(Sys.time())   
      object <- mn
      ds <- fot(mn)
      mn <- methylated(object) 
      betas(object) <- daten2 (
         mn , 
         un = unmethylated(object),
         onetwo = object@featureData@data[,ds],
         fudge, ...
     )
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
         if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
            stop('can\'t load methylumi package')
         }
         history.submitted <- as.character(Sys.time())   
         object <- mns
         ds <- fot(mns)
         betas(object) <- nasen (
            mns = methylated(object), 
            uns = unmethylated(object),
            onetwo = object@featureData@data[,ds],
            fudge
        )
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
         if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
            stop('can\'t load methylumi package')
         }
         history.submitted <- as.character(Sys.time())   
         object <- mns
         ds <- fot(mns)
         betas(object) <- dasen (
            mns = methylated(object), 
            uns = unmethylated(object),
            onetwo=mns@featureData@data[,ds],
            fudge, roco
         ) 
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
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
   history.submitted <- as.character(Sys.time())   
         object <- mns
         ds <- fot(mns)
         betas(object) <- danen (
         mns = methylated(object), 
         uns = unmethylated(object),
         onetwo=mns@featureData@data[,ds],
         fudge, ... 
       ) 
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

#tost <- function( mn, un, da, pn ) {
setMethod(
   f= "tost",
   signature(mn="MethyLumiSet"),
   definition=function(mn){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      history.submitted <- as.character(Sys.time())   
         object <- mn
         mn <- methylated(object) 
         betas(object) <- tost (
            mn, 
            un = unmethylated(object),
            da = object@featureData@data,
            pn = pvals(object)
         ) 
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
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
   history.submitted <- as.character(Sys.time())   
         object <- data
         betas(object) <- fuks (
         data = betas(object), 
         anno = object@featureData@data
       ) 
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

#swan <- function (mn, un, qc ) {
setMethod(
   f= "swan",
   signature(mn="MethyLumiSet"),
   definition=function(mn, da=NULL){
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      history.submitted <- as.character(Sys.time())   
      object <- mn
      mn <- methylated(object)
      methylumi:::betas(object) <- swan (
         mn,
         un = unmethylated(object), 
         qc = intensitiesByChannel(QCdata(object)),
         da
       ) 
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
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
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
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
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
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
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
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
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
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
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
      
      if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load methylumi package')
      }
      object <- mn
      bn     <- betas(object)
      bc     <- betas(object)
      mn     <- methylated(object)
      un     <- unmethylated(object)
      pn     <- pvals(object)
      da     <- object@featureData@data
      l      <- pfilter (
         mn=mn, un=un, bn=bn, da=da, 
         pn=pn, bc=bn, perCount, pnthresh,perc,
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
         if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
            stop('can\'t load methylumi package')
         }
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



