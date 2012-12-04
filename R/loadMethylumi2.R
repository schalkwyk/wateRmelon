loadMethylumi2 <-
function(methylationData, controlData=NULL){

	methLumi_data <- lumiMethyR2(methylationData)
	
	if(!is.null(controlData)) methLumi_data <- addControlData2methyLumiM(controlData, methLumi_data)

	return(methLumi_data)
}
