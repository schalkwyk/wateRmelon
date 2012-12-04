nbBeadsFilter <-
function(methyLumi, nbBeads.threshold = 3){

	#get nBeads for A and B signals (unmethylated and methylated signals)
	nbBeadA <- assayDataElement(methyLumi, "Avg_NBEADS_A")
	nbBeadB <- assayDataElement(methyLumi, "Avg_NBEADS_B")
	
	#set nBeadsA and B to NULL
	assayDataElement(methyLumi, "Avg_NBEADS_A") <- NULL
	assayDataElement(methyLumi, "Avg_NBEADS_B") <- NULL
	
	#identify nbBeads < 3
	indexA <- which(nbBeadA < 3)
	indexB <- which(nbBeadB < 3)
	indexAB <- union(indexA, indexB)
	rm(indexA, indexB)

	#get detection pvalues
	detectPval <- assayDataElement(methyLumi, "detection")

	#set detection p-values associated to non significant signals to 1
	detectPval[indexAB] <- 1
	assayDataElement(methyLumi, "detection") <- detectPval
	rm(detectPval, indexAB)

	return(methyLumi)
}
