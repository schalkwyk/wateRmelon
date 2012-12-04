normalizeIlluminaMethylation <-
function(
	beta,
	detect.pval,
	quantile.norm.pvalThreshold = 0.01,
	probeAnnotations,
	probeAnnotationsCategory = "relationToCpG"
	)
{
	
	#get probe IDs for each kind of Illumina bead type (Infinium I & Infinium II
	infiniumI <- probeAnnotations$TargetID[which(probeAnnotations$INFINIUM_DESIGN_TYPE == "I")]
	infiniumII <- probeAnnotations$TargetID[which(probeAnnotations$INFINIUM_DESIGN_TYPE =="II")]
	
	cat("\t Quantile normalization of samples: separated and 'robust' quantile normalization for Infinium probes I and II through probe categories (reference quantiles computed from filtered Infinium I probes only and for different categories of probe annotations). \n")
	if(!is.null(probeAnnotationsCategory)){
		if(probeAnnotationsCategory == "relationToCpG") {
			index <- which(is.element(colnames(probeAnnotations), c("TargetID", "RELATION_TO_UCSC_CPG_ISLAND")))
			probeAnnotations <- probeAnnotations[,index]
		}
		if(probeAnnotationsCategory == "relationToSequence"){
			index <- which(is.element(colnames(probeAnnotations), c("TargetID", "UCSC_REFGENE_GROUP")))
			probeAnnotations <- probeAnnotations[,index]
		}
		if(!is.element(probeAnnotationsCategory, c("relationToCpG", "relationToSequence"))){
			print("WARNINGS: probe annotation category must be one of 'relationToCpG' or 'relationToSequence'.")
			return("WARNINGS: probe annotation category must be on of 'relationToCpG' or 'relationToSequence'.")
		}
	}
	else{
		print("WARNING ! You have to specify an annotation type for probe categories based normalization ('relationToCpG' or 'relationToSequence').")
		return("WARNING ! You have to specify an annotation type for probe categories based normalization ('relationToCpG' or 'relationToSequence').")
	}

	#start subset quantile normalization, this function returns a list of 2 matrices (beta values & detection p-values)
	data.norm <- robustQuantileNorm_Illumina450K(data = beta, infiniumI = infiniumI, infiniumII = infiniumII, detect.pval = detect.pval, detect.pval.threshold=quantile.norm.pvalThreshold, annotations = probeAnnotations)
	
	names(data.norm) <- c("beta", "detection.pvalue")
	rm(beta, detect.pval)
	rm(infiniumI, infiniumII, probeAnnotations)
	
	cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1], "x", dim(data.norm$beta)[2],"\n")
	cat("Dimension of normalized detection p-values matrix: ", dim(data.norm$detection.pvalue)[1], "x", dim(data.norm$detection.pvalue)[2],"\n")
		
	return(data.norm)
}
