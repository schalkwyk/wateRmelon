#
#dfort <- function(rn){
#  # object names in IlluminaHumanMethylation450k
#  cols  <- c( "COLORCHANNEL", "CPGIRELATION", "DESIGN" )    
#  thing <- pop( cols, rn )
#  # col names in FinalReport
#  colnames(thing) <- c( "COLOR_CHANNEL", "RELATION_TO_UCSC_CPG_ISLAND", "INFINIUM_DESIGN_TYPE" )
#  data.frame( TargetID=rownames(thing), thing )
#}

dfort <- function(rn){
#  # object names in IlluminaHumanMethylation450k
#  cols  <- c( "COLORCHANNEL", "CPGIRELATION", "DESIGN" )    
stopifnot(
   all.equal(
      rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest), 
      rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Islands.UCSC)
   )
)

data.frame( 
      TargetID =   rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest),
      COLOR_CHANNEL = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest$Color,
      RELATION_TO_UCSC_CPG_ISLAND = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Islands.UCSC$Relation_to_Island,
      INFINIUM_DESIGN_TYPE =  IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest$Type
)

}


tost <-
function( mn, un, da, pn ) {
## da requirements should be checked: color channel required
# @#$%^&* mapping Illumina heads to IlluminaHumanMethylation450k.db
# 'COLOR_CHANNEL'  "IlluminaHumanMethylation450kCOLORCHANNEL"
#, 'CHROMOSOME', 'POSITION', 
# 'TargetID'  


# make a methylumi object
s <- colnames(mn)
pData <- data.frame(sampleID = s, label = s)
rownames(pData) <- s
varMetadata <- data.frame(labelDescription = colnames(pData))
rownames(varMetadata) <- colnames(pData)
data <- new("AnnotatedDataFrame", data = pData, varMetadata = varMetadata)


d               <- new("MethyLumiSet")
methylated(d)   <- mn
unmethylated(d) <- un
fData(d)        <- da
#phenoData(d)    <- as(fac, "AnnotatedDataFrame")
phenoData(d)    <- data 
pvals(d)        <- pn


e <- preprocessIlluminaMethylation(

        d,
#       path2data,
#       path2controlData,
#       projectName,
        nbBeads.threshold=NULL,
        detectionPval.threshold=NULL,
        detectionPval.perc.threshold=80,
        sample2keep =NULL,
        probeSNP_LIST=NULL,
        XY.filtering=FALSE,
        colorBias.corr=TRUE,
        bg.adjust="separatecolors",
        PATH="./"
        )

data.preprocess.norm <- normalizeIlluminaMethylation(
        beta = getMethylumiBeta(e),
        detect.pval = pvals(e),
        quantile.norm.pvalThreshold = .01,
        probeAnnotations = fData(e),
        probeAnnotationsCategory = "relationToCpG"
        )


data.preprocess.norm$beta 

}
