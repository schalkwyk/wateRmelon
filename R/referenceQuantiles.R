referenceQuantiles <-
function(data){	
	
	sortedData <- apply(data, 2, sort, na.last=TRUE)
		
	Reference.Quantiles <- rowMeans(sortedData, na.rm=TRUE)
	
	indexNA <- which(is.na(Reference.Quantiles))
	if(length(indexNA)>0) Reference.Quantiles <- Reference.Quantiles[-indexNA]
	
	return(Reference.Quantiles)
}
