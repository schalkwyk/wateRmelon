preprocessIlluminaMethylation <-
function(
        methLumi_data,
#LS#	path2data,
#LS#	path2controlData,
#LS#	projectName,
	nbBeads.threshold=NULL,
	detectionPval.threshold=NULL,
	detectionPval.perc.threshold=80,
	sample2keep,
	probeSNP_LIST,
	XY.filtering,
	colorBias.corr=TRUE,
	bg.adjust="separatecolors",
	PATH="./")
{

#set pipeline steps counter
	i=0

#loadMethylumi and store nbBeads_A and nbBeads_B info
	#LS#cat("\n\tStart data loading...\n")
	#LS#methLumi_data <- loadMethylumi2(methylationData = path2data, controlData = path2controlData)
	#LS#cat("\tProject sample nb: ", length(sampleNames(methLumi_data)), ".\n")
	#LS#cat("\tData dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
	#LS#cat("\t...data loaded..\n\n")

# nbBeads filtering
	if(!is.null(nbBeads.threshold)){
		i<-i+1
		cat(" Step ", i, ": start nb beads/probe filtering...\n")
		methLumi_data <- nbBeadsFilter(methLumi_data, nbBeads.threshold)
		cat("\t...done.\n\n")
	}
	
# sample QC an filtering
	if(!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold)){
		i<-i+1
		cat(" Step ", i, ": start sample QC & filtering...\n")
		methLumi_data <- detectionPval.filter(methLumi_data, detectionPval.threshold, detectionPval.perc.threshold, projectName, PATH=PATH)
		cat("\t Project samples nb. after sample QC: ", length(sampleNames(methLumi_data)), ".\n")
		cat("\t...done.\n\n")
	}
	
# remove controls or unrelevant samples
	if(!is.null(sample2keep)) {
		i<-i+1
		cat(" Step ", i, ": 'pertinent' sample selection...\n")
		sample2keep <- read.table(file=sample2keep, sep="\t", header=FALSE, quote="")[[1]]
		methLumi_data <- getSamples(methLumi_data, sample2keep)
		cat("\t Project samples nb after sample selection: ", length(sampleNames(methLumi_data)), ".\n")
		cat("\t...done.\n\n")
	}

# SNP filtering
	if(!is.null(probeSNP_LIST)){
		i <- i+1
		cat(" Step ", i, ": start frequent SNP filtering...\n")
		indexProbe2remove <- which(is.element(featureNames(methLumi_data), probeSNP_LIST))
		if(length(indexProbe2remove)>0) methLumi_data <- methLumi_data[-indexProbe2remove,]

		cat("\t Data dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
		cat("\t...done.\n\n")
	}
	
# XY chz filtering
	if(XY.filtering){
		i <- i+1
		cat(" Step ", i, ": start elimination of X & Y chr. probes...\n")
		methLumi_data <- filterXY(methLumi_data)
		cat("\t Data dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
		cat("\t...done.\n\n")
	}

#Color bias correction
	if(colorBias.corr){
		i <- i+1
		cat(" Step ", i, ": start color bias correction X...\n")
		methLumi_data <- lumi::adjColorBias.quantile(methLumi_data)
		cat("\t...done.\n\n")
	}

#LS#

#  I can't immediately see why the original script doesn't also require this coersion.
#  Calculates slightly different betas, approximately m/m+u ie fudge=0

 methLumi_data <-  as(methLumi_data,'MethyLumiM')
#LS#
	
#BG subtraction
	if(bg.adjust=="separatecolors"){
		i <- i+1
		cat(" Step ", i, ": start background subtraction (" , bg.adjust, ") ...\n")
		methLumi_data <- lumiMethyB(methLumi_data, separateColor = TRUE)
		cat("\t...done.\n\n")
	}
	if(bg.adjust=="unseparatecolors"){
		i <- i+1
		cat(" Step ", i, ": start background subtraction (" , bg.adjust, ") ...\n")
		methLumi_data <- lumiMethyB(methLumi_data, separateColor = FALSE)
		cat("\t...done.\n\n")
	}
	if(bg.adjust=="no"){
		i <- i+1
		cat(" Step ", i, ": no background subtraction.\n")
	}
			
	return(methLumi_data)
}
