robustQuantileNorm_Illumina450K <-
function(
	data,
	infiniumI,
	infiniumII,
	detect.pval,
	detect.pval.threshold,
	annotations,
	verbose = TRUE){
	
	sampleNames <- colnames(data)
		
	indexInfiniumI <- which(is.element(rownames(data), infiniumI))
	indexInfiniumII <- which(is.element(rownames(data), infiniumII))
	data.infiniumI <- data[indexInfiniumI,]
	data.infiniumII <- data[indexInfiniumII,]
	probeID.infiniumI <- rownames(data.infiniumI)
	probeID.infiniumII <- rownames(data.infiniumII)
	rm(data)
				
	#print(paste("Separate, robust and filtered quantile normalization per probe categories: reference quantiles for Infinium I and II probes are based on Infinium I signals associated to p-values < ", detect.pval.threshold ,". Reference quantiles are computed separately for each kind of probe annotation.", sep=""))

	data.norm <- robustQuantileNorm_Illumina450K.probeCategories(data.infiniumI, data.infiniumII, detect.pval.infiniumI = detect.pval[indexInfiniumI,], annotations = annotations, threshold = detect.pval.threshold)

	rm(data.infiniumI, data.infiniumII, indexInfiniumI, indexInfiniumII)
		
	if(verbose) {
		cat("\tNormalized beta values example:\n")
		cat("\t", data.norm[1:2,1:2],"\n")
	}
	return(coRankedMatrices(data.norm, detect.pval))
}
