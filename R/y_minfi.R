# methods for MethylSet (minfi)

#setClass("MethylSet",
# representation("VIRTUAL",
# description = "character"))

#setClass("RGChannelSet",
# representation("VIRTUAL",
# description = "character"))


# betaqn <- function (bn){ 
setMethod(
   f= "betaqn",
   signature(bn="MethylSet"),
   definition=function(bn){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object <- bn
      betaqn (
         bn = getBeta(object) 
      ) 
   }
)

setMethod(
   f= "betaqn",
   signature(bn="RGChannelSet"),
   definition=function(bn){
      if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
            stop('can\'t load minfi package')
      }
      object <- bn
      betaqn (
      bn = getBeta(object) 
      ) 
   }
)

# naten <- function (mn, un, fudge=100 ) {
setMethod(
   f= "naten",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object <- mn
      naten (
         mn = getMeth(object), 
         un = getUnmeth(object),
         fudge
      ) 
   }
)

setMethod(
   f= "naten",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mn  <- preprocessRaw(mn)    
      naten ( mn ) 
   }
)

#nanet <- function(mn, un, fudge=100){
setMethod(
   f= "nanet",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object <- mn
      nanet (
         mn = getMeth(object), 
         un = getUnmeth(object),
         fudge
      ) 
   }
)

setMethod(
   f= "nanet",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mn  <- preprocessRaw(mn)    
      nanet ( mn ) 
   }
)

got <- function(obj){
   I  <- getProbeInfo(obj)$Name
   rn <- rownames( obj@featureData@data )
   ot <- rn %in% I 
   ot[ot] <- 'I'
   ot[ot=='FALSE'] <- 'II'
   ot
}

#nanes <- function(mns, uns, onetwo, fudge=100){
setMethod(
   f= "nanes",
   signature(mns="MethylSet"),
   definition=function(mns, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
         object <- mns
         nanes (
            mns    = getMeth(object), 
            uns    = getUnmeth(object),
            onetwo = got(object),
            fudge
       ) 
   }
)

setMethod(
   f= "nanes",
   signature(mn="RGChannelSet"),
   definition=function(mns, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mns  <- preprocessRaw(mns)    
      nanes ( mns ) 
   }
)

#danes <- function(mn, un, onetwo, fudge=100, ...){ 
setMethod(
   f= "danes",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
         object <- mn
         mn <- getMeth(object)
         danes (
            mn     , 
            un     = getUnmeth(object),
            onetwo = got(object),
            fudge
       ) 
   }
)

setMethod(
   f= "danes",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mn  <- preprocessRaw(mn)    
      danes ( mn ) 
   }
)

#danet <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "danet",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
         object <- mn
         mn <- getMeth(object) 
         danet (
            mn , 
            un = getUnmeth(object),
            onetwo = got(object),
            fudge
       ) 
   }
)

setMethod(
   f= "danet",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mn  <- preprocessRaw(mn)    
      danet ( mn ) 
   }
)


#daten1 <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "daten1",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100, ...){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
         object <- mn
         mn <- getMeth(object) 
         daten1 (
            mn , 
            un = getUnmeth(object),
            onetwo = got(object),
            fudge, ...
       ) 
   }
)

setMethod(
   f= "daten1",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100, ...){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mn  <- preprocessRaw(mn)    
      daten1 ( mn, ... ) 
   }
)

#daten2 <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "daten2",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100, ...){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
         object <- mn
         mn <- getMeth(object) 
         daten2 (
            mn , 
            un = getUnmeth(object),
            onetwo = got(object),
            fudge, ...
       ) 
   }
)

setMethod(
   f= "daten2",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mn  <- preprocessRaw(mn)    
      daten2 ( mn ) 
   }
)


#ot <- function ( mns, uns, onetwo, fudge=100) {
setMethod(
   f= "nasen",
   signature(mns="MethylSet"),
   definition=function(mns, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
         object <- mns
         nasen (
            mns = getMeth(object), 
            uns = getUnmeth(object),
            onetwo = got(object),
            fudge
       ) 
   }
)

setMethod(
   f= "nasen",
   signature(mns="RGChannelSet"),
   definition=function(mns, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mns  <- preprocessRaw(mns)    
      nasen ( mns ) 
   }
)


#dasen <- function(mns, uns, onetwo, fudge=100, ...){
setMethod(
   f= "dasen",
   signature(mns="MethylSet"),
   definition=function(mns, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object <- mns
      dasen (
         mns = getMeth(object), 
         uns = getUnmeth(object),
         onetwo=got(object),
         fudge
       ) 
   }
)


setMethod(
   f= "dasen",
   signature(mns="RGChannelSet"),
   definition=function(mns, fudge=100){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mns  <- preprocessRaw(mns)    
      dasen ( mns ) 
   }
)

#danen <- function ( mns, uns, onetwo, fudge=100, ...){
setMethod(
   f= "danen",
   signature(mns="MethylSet"),
   definition=function(mns, fudge=100, ...){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
         object <- mns
         danen (
         mns = getMeth(object), 
         uns = getUnmeth(object),
         onetwo=got(object),
         fudge, ...
       ) 
   }
)

setMethod(
   f= "danen",
   signature(mns="RGChannelSet"),
   definition=function(mns, fudge=100, ...){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      mns  <- preprocessRaw(mns)    
      danen ( mns, ... ) 
   }
)


#tost <- function( mn, un, da, pn ) {  no methylset method because needs detection P values
setMethod(
   f= "tost",
   signature(mn="RGChannelSet"),
   definition=function(mn){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      pn <- detectionP(mn)
      object <- preprocessRaw(mn)
      mn <- getMeth(object) 
      tost (
         mn, 
         un = getUnmeth(object),
 #        da = object@featureData@data,  # prob need more columns
	 da = dfort(featureNames(object)),
         pn
      ) 
   }
)

#peak.correction <- function (data, anno) {
setMethod(
   f= "fuks",
   signature(data ="MethylSet"),
   definition=function(data){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object <- data
      fuks (
         data = getBeta(object), 
         anno = data.frame(DESIGN=got(object))
      ) 
   }
)

setMethod(
   f= "fuks",
   signature(data ="RGChannelSet"),
   definition=function(data){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object <- preprocessRaw(data)
      fuks (
         data = getBeta(object), 
         anno = data.frame(DESIGN=got(object))
      ) 
   }
)

#swan <- function (mn, un, qc ) {
# swan is already in minfi
setMethod(
   f= "swan",
   signature(mn ="RGChannelSet"),
   definition=function(mn){

   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }

   preprocessSWAN(mn)

})

#genki <- function(bn, g=getsnp(rownames(bn)), se=TRUE ){
setMethod(
   f= "genki",
   signature(bn="MethylSet"),
   definition=function(bn, se=TRUE){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object <- bn
      bn     <- getBeta(object)
      g      <- getsnp(rownames(bn))
      genki( bn, g, se ) 
     
   }
)

setMethod(
   f= "genki",
   signature(bn= "RGChannelSet"),
   definition=function(bn, se=TRUE){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object <- bn
      bn     <- getBeta(object)
      g      <- getsnp(rownames(bn))
      genki( bn, g, se ) 
     
   }
)


#dmrse <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse",
   signature(betas="MethylSet"),
   definition=function(betas, idmr=iDMR()){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      bn     <- getBeta(betas)
      dmrse( bn, idmr ) 
     
   }
)

setMethod(
   f= "dmrse",
   signature(betas= "RGChannelSet"),
   definition=function(betas, idmr=iDMR()){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      bn     <- getBeta(betas)
      dmrse( bn, idmr ) 
     
   }
)



#dmrse_row <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse_row",
   signature(betas="MethylSet"),
   definition=function(betas, idmr=iDMR()){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      bn     <- getBeta(betas)
      dmrse_row( bn, idmr ) 
     
   }
)

setMethod(
   f= "dmrse_row",
   signature(betas="RGChannelSet"),
   definition=function(betas, idmr=iDMR()){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      bn     <- getBeta(betas)
      dmrse_row( bn, idmr ) 
     
   }
)


#dmrse_col <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse_col",
   signature(betas="MethylSet"),
   definition=function(betas, idmr=iDMR()){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      bn     <- getBeta(betas)
      dmrse_col( bn, idmr ) 
     
   }
)

setMethod(
   f= "dmrse_col",
   signature(betas="RGChannelSet"),
   definition=function(betas, idmr=iDMR()){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      bn     <- getBeta(betas)
      dmrse_col( bn, idmr ) 
     
   }
)

#seabi <- function (bn, stop, sex, X){
setMethod(
   f= "seabi",
   signature(bn="MethylSet"),
   definition=function( bn, stop=1, sex, X ){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      betas    <- getBeta(bn)
      seabi( betas, stop, sex, X ) 
     
   }
)

setMethod(
   f= "seabi",
   signature(bn="RGChannelSet"),
   definition=function( bn, stop=1, sex, X ){
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      betas    <- getBeta(bn)
      seabi( betas, stop, sex, X ) 
     
   }
)

#pfilter<-function(mn, un, bn, da, onetwo, pn, bc, perCount, pnthresh, perc, pthresh){    
# filter function by Ruth Pidsley

setMethod(
   f= "pfilter",
   signature(mn="RGChannelSetExtended"),
   definition=function(
      mn=RGChannelSetExtended, 
      perCount =NULL, pnthresh=NULL, perc=NULL, 
      pthresh=NULL
   ){
      if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
      }
      object   <- mn
      object2  <- preprocessRaw(object)
      mn       <- getMeth(object2)
      un       <- getUnmeth(object2)
      pn       <- detectionP(object)
      bc       <- beadcount(object)
      l        <- pfilter (
         mn=mn, un=un, bn=NULL, 
         da=NULL, 
         pn=pn, bc=bc,perCount, 
	 pnthresh,perc,pthresh
      ) 
      object3 <-new("MethylSet",
         Meth=l$mn,
         Unmeth=l$un,
         annotation=annotation(object2),
         phenoData=phenoData(object2[,colnames(l$mn)])
      )  			      	
      object3@preprocessMethod <- c(
         "Raw (no normalization or bg correction) with wateRmelon pfilter",
         as.character(packageVersion("minfi")), 
         as.character(packageVersion("IlluminaHumanMethylation450kmanifest"))
      )
      object3     
   }
)


