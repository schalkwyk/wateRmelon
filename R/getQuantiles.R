getQuantiles <-
function(source, targetLength){

	if(targetLength < 2){
		print("WARNING! in getQuantiles.R, 'targetLenght' must be > 1")
		return("WARNING! in getQuantiles.R, 'targetLenght' must be > 1")
	}

	p <- 1/(targetLength-1)

	probs <- vector(, targetLength)
	probs[1]=0

	for(i in 2:targetLength) probs[i] <- probs[i-1]+p
	
	probs[targetLength] <- 1

	targetQuantiles <- quantile(x=sort(source), probs=probs, na.rm=TRUE)

	return(targetQuantiles)
}
