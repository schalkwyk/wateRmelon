getSamples <-
function(methLumi_data, sample2keep){

	samples <- sampleNames(methLumi_data)
	index <- which(is.element(samples, sample2keep))

	if(length(index)>0) methLumi_data <- methLumi_data[, index]
	else(print("WARNING ! No samples to select !"))
	
	return(methLumi_data)
}
