lumiMethyR2 <-
function (data, lib = NULL, controlData = NULL){

    methyLumiSet <- methylumiR(data)

	#sort by sample names
	indexSamples <- sort(sampleNames(methyLumiSet), index.return=TRUE)$ix
	methyLumiSet <- methyLumiSet[,indexSamples]
	
	methyLumiM <- as(methyLumiSet, "MethyLumiM")
	# 'A' for unmethylated and 'B' for methylated
	assayDataElement(methyLumiM, "Avg_NBEADS_A") <- assayDataElement(methyLumiSet, "Avg_NBEADS_A")
	assayDataElement(methyLumiM, "Avg_NBEADS_B") <- assayDataElement(methyLumiSet, "Avg_NBEADS_B")
	rm(methyLumiSet)
	
    if (!is.null(lib)) methyLumiM <- addAnnotationInfo(methyLumiM, lib = lib)
    if (!is.null(controlData)){
        if (is.character(controlData)) controlData <- methylumiR(controlData)
        if (is(controlData, "MethyLumiQC")) controlData(methyLumiM) <- controlData
        else { cat("Provided controlData is not supported!\n") }
    }
	
    return(methyLumiM)
}
