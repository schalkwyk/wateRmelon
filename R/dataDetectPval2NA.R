dataDetectPval2NA <-
function(data, detect.pval, threshold){

	indexNA <- which(detect.pval > threshold, arr.ind=TRUE)

	data[indexNA] <- NA
					
	return(data)
}
