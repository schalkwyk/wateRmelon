#' Functions from 450-pipeline (Touleimat & Tost)
#' 
#' These functions are part of the 450K pipeline (Touleimat and Tost,
#' Epigenomics 2012 4:325).  For freestanding use of the normalization
#' function, a wrapper is provided, see \code{\link{tost}}
#' 
#' 
#' @aliases adaptRefQuantiles coRankedMatrices concatenateMatrices
#' dataDetectPval2NA detectionPval.filter filterXY findAnnotationProbes
#' getMethylumiBeta getQuantiles getSamples loadMethylumi2 lumiMethyR2
#' nbBeadsFilter normalize.quantiles2 normalizeIlluminaMethylation
#' pipelineIlluminaMethylation.batch preprocessIlluminaMethylation
#' referenceQuantiles robustQuantileNorm_Illumina450K
#' robustQuantileNorm_Illumina450K.probeCategories uniqueAnnotationCategory
#' bgIntensitySwan.methylumi
#' @return see \code{\link{tost}} %% ~Describe the value returned %% If it is a
#' LIST, use %% \item{comp1 }{Description of 'comp1'} %% \item{comp2
#' }{Description of 'comp2'} %% ...
#' @author Nizar Touleimat, wrapper by lschal@@essex.ac.uk
#' @references Touleimat N, Tost J: Complete pipeline for Infinium R Human
#' Methylation 450K BeadChip data processing using subset quantile
#' normalization for accurate DNA methylation estimation. Epigenomics 2012,
#' 4:325-341
#' 
#' Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A data-driven
#' approach to preprocessing Illumina 450K methylation array data (submitted)
#' @export adaptRefQuantiles
adaptRefQuantiles <-
function(Reference.Quantiles, sizeNew, verbose=TRUE){
	#print(paste("Length sizeNew2: ", sizeNew, sep=""))

	M <- length(Reference.Quantiles)

	if(sizeNew == M){
		newReference.Quantiles <- Reference.Quantiles
		if(verbose) cat("\t\t adaptRefQuantiles:  sizeNew = length(Reference.Quantiles)")
	}
	else{ newReference.Quantiles <- getQuantiles(na.exclude(Reference.Quantiles), sizeNew)
		if(verbose){
			if(sizeNew < M) cat("\t\t adaptRefQuantiles:  sizeNew < length(Reference.Quantiles)")
			if(sizeNew > M) cat("\t\t adaptRefQuantiles:  sizeNew > length(Reference.Quantiles)")
		}
	}
	return(sort(newReference.Quantiles))
}
