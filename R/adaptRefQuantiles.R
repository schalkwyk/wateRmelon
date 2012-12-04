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
