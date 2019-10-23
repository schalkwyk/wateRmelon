pipelineIlluminaMethylation.batch <-
function(
	PATH_PROJECT_DATA,
	projectName,
	nbBeads.threshold,
	detectionPval.threshold,
	detectionPval.perc.threshold,
	sampleSelection,
	probeSNP_LIST,
	XY.filtering,
	colorBias.corr,
	bg.adjust,
	PATH
 ){

	####################################
	# Get project pathes and load data #
	####################################

	subProjects <- dir(PATH_PROJECT_DATA)

	#for all subprojects
	for(i in 1:length(subProjects)){
		projectName <- subProjects[i]
		sampleTable <- dir(paste(PATH_PROJECT_DATA, projectName, "/", sep=""), pattern="TableSample")
		controlTable <- dir(paste(PATH_PROJECT_DATA, projectName, "/", sep=""), pattern="TableControl")
		cat("\n# ")
		cat("Processing sub-project: ", projectName, "\n")

		if(length(sampleTable) > 1){
			warning <- "\tWARNING ! Sample table: too many files matching with pattern 'TableSample' ! \n"
			cat(warning)
			return(warning)
		}
		if(length(sampleTable) < 1){
			warning <- "\tWARNING ! Sample table: no file matching with pattern 'TableSample' ! \n"
			cat(warning)
			return(warning)
		}
		if(length(controlTable) > 1){
			warning <- "\tWARNING ! Control table: too many files matching with pattern 'TableControl' ! \n"
			cat(warning)
			return(warning)
		}
		if(length(controlTable) < 1){
			warning <- "\tWARNING ! Control table: no file matching with pattern 'TableControl' ! \n"
			cat(warning)
			return(warning)
		}

		if(sampleSelection){
			sampleList <- dir(paste(PATH_PROJECT_DATA, projectName, "/", sep=""), pattern = "sampleList")
			path2sampleList <- paste(PATH_PROJECT_DATA, projectName, "/", sampleList, sep="")

			if(length(sampleList) > 1){
				warning <- "\tWARNING ! List for sample selection: too many files matching with pattern 'SampleList' ! \n"
				cat(warning)
				return(warning)
			}
			if(length(sampleList) < 1){
				warning <- "\tWARNING ! List for sample selection: no file matching with pattern 'SampleList' ! \n"
				path2sampleList <- NULL
				cat(warning)
			}
		}
		else{path2sampleList <- NULL}

		#path2data <- paste(PATH_PROJECT_DATA, projectName, "/", sampleTable, sep="")
		path2data <-  sampleTable
		path2controlData <- paste(PATH_PROJECT_DATA, projectName, "/", controlTable, sep="")

		cat("\tSample table: ", path2data, "\n")
		cat("\tControl table: ", path2controlData, "\n")
		cat("\tSample list (for sample selection): ", path2sampleList, "\n")

		#############################
		# starts data preprocessing #
		#############################

		methLumi_data <- preprocessIlluminaMethylation(
		# TGS	path2data = path2data,
		# TGS	path2controlData = path2controlData,
		# TGS	projectName = projectName,
			nbBeads.threshold = nbBeads.threshold,
			detectionPval.threshold = detectionPval.threshold,
			detectionPval.perc.threshold = detectionPval.perc.threshold,
			sample2keep = path2sampleList,
			probeSNP_LIST,
			XY.filtering = XY.filtering,
			colorBias.corr = colorBias.corr,
			bg.adjust = bg.adjust,
			PATH = PATH_RES
		)

		################################################
		# Sub-project data & information concatenation #
		################################################

		if(i == 1){
			beta <- getMethylumiBeta(methLumi_data)
			cat("\t beta plate", i, " ok (", dim(beta)[1], "x", dim(beta)[2], ").\n")
			detectionPval <- assayDataElement(methLumi_data, "detection")
			cat("\t detection p-values plate", i, " ok (", dim(detectionPval)[1], "x", dim(detectionPval)[2], ").\n")
			#select "useful" probe annotations
			annotation <- fData(methLumi_data) ; rm(methLumi_data)
			index <- which(is.element(colnames(annotation), c("TargetID", "INFINIUM_DESIGN_TYPE", "RELATION_TO_UCSC_CPG_ISLAND", "UCSC_REFGENE_GROUP")))
			annotation <- annotation[,index]
		}
		#concatenate 'betas'
		else{
			beta_i <- getMethylumiBeta(methLumi_data)
			cat("\t beta_", i, " ok (", dim(beta_i)[1], "x", dim(beta_i)[2], ").\n")
			detectionPval_i <- assayDataElement(methLumi_data, "detection")
			cat("\t For all sub-projects: beta matrices concatenation & detection p-value matrices concatenation.\n")

			beta <- concatenateMatrices(beta, beta_i) ; rm(beta_i)
			detectionPval <- concatenateMatrices(detectionPval, detectionPval_i) ; rm(detectionPval_i)
			annotation <- annotation[ which(is.element(annotation$TargetID, rownames(beta))),]

			cat("\t beta ok (", dim(beta)[1], "x", dim(beta)[2], ").\n")
			cat("\t detection p-values ok (", dim(detectionPval)[1], "x", dim(detectionPval)[2], ").\n")
		}
	}

	############################################################################################
	# Extraction of SNP probes ("rs" probes)
	############################################################################################

	indexSNP <- grep(pattern="rs*", x=rownames(beta))
	if(length(indexSNP)>0){
		betaSNP <- beta[indexSNP,]
		detectionPvalSNP <- detectionPval[indexSNP,]

		beta <- beta[-indexSNP,]
		detectionPval <- detectionPval[-indexSNP,]

		write.csv(betaSNP, file=paste(PATH_RES, projectName, "_betaSNPprobes.csv", sep=""), quote=FALSE)
		write.csv(detectionPvalSNP, file=paste(PATH_RES, projectName, "_detectionPvalueSNPprobes.csv", sep=""), quote=FALSE)
	}

	############################################################################################
	# start data normalization (subset quantile normalization per probe annotation categories) #
	############################################################################################

	data.preprocess.norm <- normalizeIlluminaMethylation(
		beta = beta,
		detect.pval = detectionPval,
		quantile.norm.pvalThreshold = detectionPval.threshold,
		probeAnnotations = annotation,
		probeAnnotationsCategory = probeAnnotationsCategory
	)

	return(data.preprocess.norm)
}
