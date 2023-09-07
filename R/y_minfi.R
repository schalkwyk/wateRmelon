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
      object <- mn
      out <- naten (
         mn = getMeth(object),
         un = getUnmeth(object),
         fudge, ret2=T
      )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(mn),
                     annotation = annotation(mn), metadata = metadata(mn))
   out2@preprocessMethod <- c(rg.norm = 'naten (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(mn@annotation))))
   return(out2)
   }
)

setMethod(
   f= "naten",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
      mn  <- preprocessRaw(mn)
      naten ( mn )
   }
)

#nanet <- function(mn, un, fudge=100){
setMethod(
   f= "nanet",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100){
      object <- mn
      out <- nanet (
         mn = getMeth(object),
         un = getUnmeth(object),
         fudge, ret2=T
      )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(object),
                     annotation = annotation(object), metadata = metadata(object))
   out2@preprocessMethod <- c(rg.norm = 'nanet (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(object@annotation))))
   return(out2)
   }
)

setMethod(
   f= "nanet",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
      mn  <- preprocessRaw(mn)
      nanet ( mn )
   }
)



#' Internal functions for Illumina i450 normalization functions
#'
#' got and fot find the annotation column differentiating type I and type II
#' assays in MethylSet (got) or MethyLumiSet (fot) objects. pop extracts
#' columns from IlluminaHumanMethylation450k.db
#'
#' \code{got} returns a character vector of 'I' and 'II', \code{fot} returns
#' the index of the relevant column. \code{pop} returns a data frame
#'
#' @aliases got fot pop
#' @param x a MethyLumiSet
#' @param obj a MethylSet
#' @param fd a character vector of the desired annotation columns
#' @param rn a character vector of the desired features
#' @author lschal@@essex.ac.uk
#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted)
#' @export got

# 'Borrow Internal, unexported, function from minfi'
.getManifestString <- function(annotation) {
    if(length(annotation) == 1)
        return(paste0(annotation, "manifest"))
    if("array" %in% names(annotation))
        return(paste0(annotation["array"], "manifest"))
    stop("unable to get the manifest string for this object")
}

got <- function(obj){
   I  <- getProbeInfo(obj)$Name
   rn <- rownames(obj) # FeatureData is removed from MethylSet
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
         object <- mns
         out <- nanes (
            mns    = getMeth(object),
            uns    = getUnmeth(object),
            onetwo = got(object),
            fudge,
            ret2   = TRUE
       )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(mns),
                     annotation = annotation(mns), metadata = metadata(mns))
   out2@preprocessMethod <- c(rg.norm = 'nanes (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(mns@annotation))))
   return(out2)
   }
)

setMethod(
   f= "nanes",
   signature(mn="RGChannelSet"),
   definition=function(mns, fudge=100){
      mns  <- preprocessRaw(mns)
      nanes ( mns )
   }
)

#danes <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "danes",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100){
         object <- mn
         mn <- getMeth(object)
         out <- danes (
            mn     ,
            un     = getUnmeth(object),
            onetwo = got(object),
            fudge, ret2=T
       )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(object),
                     annotation = annotation(object), metadata = metadata(object))
   out2@preprocessMethod <- c(rg.norm = 'danes (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(object@annotation))))
   return(out2)
   }
)

setMethod(
   f= "danes",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
      mn  <- preprocessRaw(mn)
      danes ( mn )
   }
)

#danet <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "danet",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100){
         object <- mn
         mn <- getMeth(object)
         out <- danet (
            mn ,
            un = getUnmeth(object),
            onetwo = got(object),
            fudge, ret2=T
       )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(object),
                     annotation = annotation(object), metadata = metadata(object))
   out2@preprocessMethod <- c(rg.norm = 'danet (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(object@annotation))))
   return(out2)
   }
)

setMethod(
   f= "danet",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
      mn  <- preprocessRaw(mn)
      danet ( mn )
   }
)


#daten1 <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "daten1",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100, ...){
         object <- mn
         mn <- getMeth(object)
         out <- daten1 (
            mn ,
            un = getUnmeth(object),
            onetwo = got(object),
            fudge, ret2=T, ...
       )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(object),
                     annotation = annotation(object), metadata = metadata(object))
   out2@preprocessMethod <- c(rg.norm = 'daten1 (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(object@annotation))))
   return(out2)
   }
)

setMethod(
   f= "daten1",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100, ...){
      mn  <- preprocessRaw(mn)
      daten1 ( mn, ... )
   }
)

#daten2 <- function(mn, un, onetwo, fudge=100, ...){
setMethod(
   f= "daten2",
   signature(mn="MethylSet"),
   definition=function(mn, fudge=100, ...){
         object <- mn
         mn <- getMeth(object)
         out <- daten2 (
            mn ,
            un = getUnmeth(object),
            onetwo = got(object),
            fudge, ret2=T, ...
       )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(object),
                     annotation = annotation(object), metadata = metadata(object))
   out2@preprocessMethod <- c(rg.norm = 'daten2 (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(object@annotation))))
   return(out2)
   }
)

setMethod(
   f= "daten2",
   signature(mn="RGChannelSet"),
   definition=function(mn, fudge=100){
      mn  <- preprocessRaw(mn)
      daten2 ( mn )
   }
)


#ot <- function ( mns, uns, onetwo, fudge=100) {
setMethod(
   f= "nasen",
   signature(mns="MethylSet"),
   definition=function(mns, fudge=100){
         object <- mns
         out <- nasen (
            mns = getMeth(object),
            uns = getUnmeth(object),
            onetwo = got(object),
            fudge, ret2=T
       )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(mns),
                     annotation = annotation(mns), metadata = metadata(mns))
   out2@preprocessMethod <- c(rg.norm = 'nasen (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(mns@annotation))))
   return(out2)
   }
)

setMethod(
   f= "nasen",
   signature(mns="RGChannelSet"),
   definition=function(mns, fudge=100){
      mns  <- preprocessRaw(mns)
      nasen ( mns )
   }
)


#dasen <- function(mns, uns, onetwo, fudge=100, ...){
setMethod(
   f= "dasen",
   signature(mns="MethylSet"),
   definition=function(mns, fudge=100){
      object <- mns
      out <- dasen (
         mns = getMeth(object),
         uns = getUnmeth(object),
         onetwo=got(object),
         fudge, ret2=T
       )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(mns),
                     annotation = annotation(mns), metadata = metadata(mns))
   out2@preprocessMethod <- c(rg.norm = 'dasen (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(mns@annotation))))
   return(out2)
   }
)


setMethod(
   f= "dasen",
   signature(mns="RGChannelSet"),
   definition=function(mns, fudge=100){
      mns  <- preprocessRaw(mns)
      dasen ( mns )
   }
)

#danen <- function ( mns, uns, onetwo, fudge=100, ...){
setMethod(
   f= "danen",
   signature(mns="MethylSet"),
   definition=function(mns, fudge=100, ...){
         object <- mns
         out <- danen (
         mns = getMeth(object),
         uns = getUnmeth(object),
         onetwo=got(object),
         fudge, ret2=T,...
       )
   out2 <- MethylSet(Meth=out[[1]], Unmeth = out[[2]], colData = colData(object),
                     annotation = annotation(object), metadata = metadata(object))
   out2@preprocessMethod <- c(rg.norm = 'danen (wateRmelon)',
                              minfi = as.character(packageVersion('minfi')),
                              manifest = as.character(packageVersion(.getManifestString(object@annotation))))
   return(out2)
   }
)

setMethod(
   f= "danen",
   signature(mns="RGChannelSet"),
   definition=function(mns, fudge=100, ...){
      mns  <- preprocessRaw(mns)
      danen ( mns, ... )
   }
)


setMethod(
   f= "uSexQN",
   signature(mns="RGChannelSet"),
   definition = function(mns, cores=1, fudge=100,...){
      mns  <- preprocessRaw(mns)
      return(uSexQN(mns, cores=cores, fudge=fudge, ...))
   }
)

setMethod(
   f= "uSexQN",
   signature(mns="MethylSet"),
   definition = function(mns, cores=1, fudge=100,...){
      object <- mns
      out <- uSexQN(
         mns = getMeth(object),
         uns = getUnmeth(object),
         ot = got(object),
         chr = getAnnotation(object)[rownames(object),'chr'],
         cores=cores,
         fudge=fudge,
         ret2=TRUE
      )
      out2 <- MethylSet(
         Meth = out$methylated,
         Unmeth = out$unmethylated,
         colData = colData(object),
         annotation = annotation(object),
         metadata = metadata(object)
      )
      out2@preprocessMethod <- c(rg.norm = 'uSexQN (wateRmelon)',
                                 minfi = as.character(packageVersion('minfi')),
                                 manifest = as.character(packageVersion(.getManifestString(object@annotation))))
      return(out2)
   }
)







#tost <- function( mn, un, da, pn ) {  no methylset method because needs detection P values
setMethod(
   f= "tost",
   signature(mn="RGChannelSet"),
   definition=function(mn){
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


   preprocessSWAN(mn)

})

#genki <- function(bn, g=getsnp(rownames(bn)), se=TRUE ){
setMethod(
   f= "genki",
   signature(bn="MethylSet"),
   definition=function(bn, se=TRUE){
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
#     object <- bn
#     bn     <- getBeta(object)
#     g      <- getsnp(rownames(bn))
#     genki( bn, g, se )

      object <- bn
      bn     <- getSnpBeta(object)
      g      <- getsnp(rownames(bn))
      genki( bn, g, se )
   }
)


#dmrse <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse",
   signature(betas="MethylSet"),
   definition=function(betas, idmr=iDMR()){
      bn     <- getBeta(betas)
      if (length(grep('_', head(rownames(bn), n = 10L)))==10){
          bn <- epicv2clean(bn)
      }
      
      dmrse( bn, idmr )

   }
)

setMethod(
   f= "dmrse",
   signature(betas= "RGChannelSet"),
   definition=function(betas, idmr=iDMR()){
      bn     <- getBeta(betas)
      if (length(grep('_', head(rownames(bn), n = 10L)))==10){
          bn <- epicv2clean(bn)
      }
      dmrse( bn, idmr )

   }
)



#dmrse_row <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse_row",
   signature(betas="MethylSet"),
   definition=function(betas, idmr=iDMR()){
      bn     <- getBeta(betas)
      if (length(grep('_', head(rownames(bn), n = 10L)))==10){
          bn <- epicv2clean(bn)
      }
      dmrse_row( bn, idmr )

   }
)

setMethod(
   f= "dmrse_row",
   signature(betas="RGChannelSet"),
   definition=function(betas, idmr=iDMR()){
      bn     <- getBeta(betas)
      if (length(grep('_', head(rownames(bn), n = 10L)))==10){
          bn <- epicv2clean(bn)
      }
      dmrse_row( bn, idmr )

   }
)


#dmrse_col <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse_col",
   signature(betas="MethylSet"),
   definition=function(betas, idmr=iDMR()){
      bn     <- getBeta(betas)
      if (length(grep('_', head(rownames(bn), n = 10L)))==10){
          bn <- epicv2clean(bn)
      }
      dmrse_col( bn, idmr )

   }
)

setMethod(
   f= "dmrse_col",
   signature(betas="RGChannelSet"),
   definition=function(betas, idmr=iDMR()){
      bn     <- getBeta(betas)
      if (length(grep('_', head(rownames(bn), n = 10L)))==10){
          bn <- epicv2clean(bn)
      }
      dmrse_col( bn, idmr )

   }
)

#seabi <- function (bn, stop, sex, X){
setMethod(
   f= "seabi",
   signature(bn="MethylSet"),
   definition=function( bn, stop=1, sex, X ){
      betas    <- getBeta(bn)
      seabi( betas, stop, sex, X )

   }
)

setMethod(
   f= "seabi",
   signature(bn="RGChannelSet"),
   definition=function( bn, stop=1, sex, X ){
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
      object   <- mn
      object2  <- preprocessRaw(object)
#      mn       <- getMeth(object2)
#      un       <- getUnmeth(object2)
      pn       <- detectionP(object)
      bc       <- beadcount(object)
      l        <- pfilter (
         mn=NULL, un=NULL, bn=NULL,
         da=NULL,
         pn=pn, bc=bc,
         perCount,
	 pnthresh,perc,pthresh,
         logical.return= TRUE
      )

      include <- names(which(l$probes))
      ret <- subsetByLoci(object, includeLoci = include)[,l$samples]
      return(ret)
    }
)

#BMIQ <- function(beta.v,design.v,nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1){
setMethod(
   f= "BMIQ",
   signature(beta.v="MethylSet"),
   definition=function(
		 beta.v,
      nL=3, doH=TRUE, nfit=5000,
      th1.v=c(0.2,0.75), th2.v=NULL,
      niter=5, tol=0.001, plots=FALSE,
      pri=FALSE

		       ){
      object <- beta.v
      d <- as.numeric(factor(got(object)))
      ibetas <- getBeta(object)
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
   }
)
# bscon <- function(x, ...) { UseMethod (bscon, x ) }
# see also bscon_methy and bscon_minfi

setMethod(
   f= "bscon",
   signature(x="RGChannelSet"),
   definition=function( x ){
      bscon_minfi(x)
   }
)

# outlyx <- function(x, y, dist1, dist2)
setMethod(
   f= "outlyx",
   signature(x="RGChannelSet"),
   definition=function(x, iqr, iqrP, pc, mv, mvP, plot){
   x <- getBeta(x)
   outlyx(x, iqr, iqrP, pc, mv, mvP, plot)
   }
)

setMethod(
   f= "outlyx",
   signature(x="MethylSet"),
   definition=function(x, iqr, iqrP, pc, mv, mvP, plot){
   x <- getBeta(x)
   outlyx(x, iqr, iqrP, pc, mv, mvP, plot)
   }
)

setMethod(
   f= "pwod",
   signature(object="RGChannelSet"),
   definition=function(object, mul){
   object <- getBeta(object)
   pwod(object, mul)
   }
)

setMethod(
   f= "pwod",
   signature(object="MethylSet"),
   definition=function(object, mul){
   object <- getBeta(object)
   pwod(object, mul)
   }
)

setMethod(
  f= "agep",
  signature(betas="MethylSet"),
  definition=function(betas, coeff, method='horvath'){
    object <- getBeta(betas)
    agep(betas=object, coeff, method=method)
  }
)

setMethod(
  f= "estimateCellCounts.wateRmelon",
  signature(object="RGChannelSet"),
  definition=function(object, referencePlatform = NULL, mn = NULL, un = NULL, bn = NULL, 
		      perc = 1, compositeCellType = "Blood", probeSelect = "auto", cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
                      returnAll = FALSE, meanPlot = FALSE, verbose = TRUE, ...){
    if(is.null(referencePlatform)){
	referencePlatform <- getManifest(object)@annotation
    }
    object <- preprocessRaw(object)
    estimateCellCounts.wmln(object=object, referencePlatform=referencePlatform, mn=getMeth(object), un=getUnmeth(object), bn=getBeta(object), 
			    perc=perc, compositeCellType=compositeCellType,
			    probeSelect=probeSelect, cellTypes=cellTypes, returnAll=returnAll, meanPlot=meanPlot, verbose=verbose)
  }
)

setMethod(
  f= "estimateCellCounts.wateRmelon",
  signature(object="MethylSet"),
  definition=function(object, referencePlatform = NULL, mn = NULL, un = NULL, bn = NULL, 
		      perc = 1, compositeCellType = "Blood", probeSelect = "auto", cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
                      returnAll = FALSE, meanPlot = FALSE, verbose = TRUE, ...){
    if(is.null(referencePlatform)){
	referencePlatform <- getManifest(object)@annotation
    }
    estimateCellCounts.wmln(object=object, referencePlatform=referencePlatform, mn=getMeth(object), un=getUnmeth(object), bn=getBeta(object), 
			    perc=perc, compositeCellType=compositeCellType, probeSelect=probeSelect, cellTypes=cellTypes, 
			    returnAll=returnAll, meanPlot=meanPlot, verbose=verbose)
  }
)



setMethod(
  f= "adjustedDasen",
  signature(mns="MethylSet"),
  definition = function(mns, offset_fit=TRUE, cores=1, fudge=100, ...){
      object <- mns
      chr <- as.character(.createAnnotation(object)$chr)
      out <- adjustedDasen(
         mns=getMeth(object),
         uns = getUnmeth(object),
         onetwo = got(object),
         chr = getAnnotation(object)[rownames(object),'chr'],
         fudge=fudge,
         cores=cores,
         offset_fit=offset_fit,
         ret2=TRUE
      )
      out2 <- MethylSet(
         Meth = out$methylated,
         Unmeth = out$unmethylated,
         colData = colData(object),
         annotation = annotation(object),
         metadata = metadata(object)
      )
      out2$preprocessMethod <- c(rg.norm = ifelse(offset_fit, 'adjustedDasen (wateRmelon)', 'adjustedNasen (wateRmelon)'),
                                 minfi = as.character(packageVersion('minfi')),
                                 manifest = as.character(packageVersion(.getManifestString(object@annotation))))
      return(out)
  }
)

setMethod(
  f= "adjustedDasen",
  signature(mns="RGChannelSet"),
  definition = function(mns, offset_fit=TRUE, cores=1, fudge=100, ...){
      object  <- preprocessRaw(mns)
      return(adjustedDasen(object, cores = cores, fudge = fudge, offset_fit = offset_fit, ...))
  }
)

