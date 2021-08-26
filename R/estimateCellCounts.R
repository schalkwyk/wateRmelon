.normalizeQuantiles2 <- function(A) {
    # Abridged function of limma::normalizeQuantiles, cuts function
    # after the sort.  Strictly used to generate quantiles.
    # A = a matrix
	n <- dim(A)
	if(is.null(n)) return(A)
	if(n[2]==1) return(A)
	O <- S <- array(,n)
	nobs <- rep(n[1],n[2])
	i <- (0:(n[1]-1))/(n[1]-1)
	for (j in 1:n[2]) {
		Si <- sort(A[,j], method="quick", index.return=TRUE)
		nobsj <- length(Si$x)
		if(nobsj < n[1]) {
			nobs[j] <- nobsj
			isna <- is.na(A[,j])
			S[,j] <- approx((0:(nobsj-1))/(nobsj-1), Si$x, i, ties="ordered")$y
			O[!isna,j] <- ((1:n[1])[!isna])[Si$ix]
		} else {
			S[,j] <- Si$x
			O[,j] <- Si$ix
		}
	}
m <- rowMeans(S)
output <- list(m, i)
return(output)
}

.impose <- function(matrix, quan){
    # Other half of normalizeQuantiles, assumes ties = TRUE
    # Used to fill in a given matrix, according to quantiles determined by quan
    # matrix = matrix such as those obtained by betas(object)
    # quan = a list of 4 vectors corresponding to quantiles, names, intervals 
    # and probe design specific to IlluminaDNAMethylationMicro-arrays.
    ot <- quan[['onetwo']]
    quantiles <- quan[['quantiles']]
    names(ot) <- quan[['rn']]
    inter <- quan[['inter']]
    blank <- matrix(NA, length(ot), ncol(matrix))
    rownames(blank) <- names(ot)
    share <- intersect(rownames(blank), rownames(matrix))
    blank[share,] <- matrix[share, ]
    for(z in 1:ncol(matrix)){
        isna <- is.na(blank[,z])
        r <- rep(0, length(ot))
        r[ot=='I'] <- rank(blank[ot=='I', z])
        r[ot=='II'] <- rank(blank[ot=='II', z])
        blank[ot == 'I' & (!isna), z] <-  approx(inter[ot == 'I'],
            quantiles[ot == 'I'],
            (r[ot=='I'&(!isna)] - 1) /(sum(ot == 'I')-1), ties = "ordered")$y

        blank[ot == 'II' & (!isna),z] <- approx(inter[ot == 'II'],
            quantiles[ot == 'II'],
            (r[ot=='II'&(!isna)] - 1)/(sum(ot == 'II')-1), ties = "ordered")$y
    }
    b2 <- na.omit(blank)
    b2 <- b2[share,]
    return(b2)
}

estimateCellCounts.wmln <- function(
    object,
    referencePlatform = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k"),
    mn = NULL,
    un = NULL,
    bn = NULL,
    perc = 1,
    compositeCellType = "Blood",
    probeSelect = "auto",
    cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),

    returnAll = FALSE,
    meanPlot = FALSE,
    verbose = TRUE,
    ...) {
    # Mostly a copy of minfi::estimateCellCounts with partial optimisation
    # to make use of some improvements that improve speed at a marginal cost of
    # accuracy. Particularly useful for big data where normalising biological
    # and reference data _together_ is as important.
    referencePlatform <- match.arg(referencePlatform)
    rgPlatform <- platform <- gsub('IlluminaHumanMethylation', '', referencePlatform)
    # Sanity Checking from minfi...
    if((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes)){
        message("[estimateCellCounts] Consider including 'nRBC' in argument 'cellTypes' for cord blood estimation.\n")
    }
    #FlowSorted.Blood.EPIC
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, rgPlatform)
    if(!require(referencePkg, character.only = TRUE)){
        stop(sprintf("Could not find reference data package for compositeCellType '%s' and referencePlatform '%s' (inferred package name is '%s')",
                     compositeCellType, rgPlatform, referencePkg))
    }
    if(platform == 'EPIC'){
	# This is actually worse than the original grok...
	if(require('ExperimentHub', character.only = TRUE)){ 
	  hub <- ExperimentHub::ExperimentHub()  
          query(hub, "FlowSorted.Blood.EPIC")
	  referenceRGset <- hub[["EH1136"]]
        } else {
	  stop('Could not find ExperimentHub R package - please install')
	}
	cellTypes <- gsub('Gran', 'Neu', cellTypes)
    } else {
    	data(list = referencePkg)
    	referenceRGset <- get(referencePkg)
    }
    if(! "CellType" %in% names(colData(referenceRGset)))
        stop(sprintf("the reference sorted dataset (in this case '%s') needs to have a phenoData column called 'CellType'"),
             names(referencePkg))
    if(sum(colnames(object) %in% colnames(referenceRGset)) > 0)
        stop("the sample/column names in the user set must not be in the reference data ")
    if(!all(cellTypes %in% referenceRGset$CellType))
        stop(sprintf("all elements of argument 'cellTypes' needs to be part of the reference phenoData columns 'CellType' (containg the following elements: '%s')",
                     paste(unique(referenceRGset$cellType), collapse = "', '")))
    if(length(unique(cellTypes)) < 2)
        stop("At least 2 cell types must be provided.")

    # Here is where things get different.
    if(is.null(mn)) mn <- methylated(object)
    if(is.null(un)) un <- unmethylated(object)
    if(is.null(bn)) bn <- betas(object)
    referencePd <- colData(referenceRGset)
    referenceMset <- preprocessRaw(referenceRGset)
    mrn <- intersect(rownames(referenceMset), rownames(mn))
    referenceMset <- referenceMset[mrn,]
    M <- mn[mrn,]
    U <- un[mrn,]
    #ot <- getProbeType(referenceMset)
    pri <- rbind(
        cbind(getProbeInfo(referenceMset, type='I')[,1], type='I'),
        cbind(getProbeInfo(referenceMset, type='II')[,1], type='II')
    )
    ot <- pri[,2]
    names(ot) <- pri[,1]
    ot <- ot[mrn]
    colsel <- sample(seq_len(ncol(M)), max(2, min(ncol(M),round(ncol(M)*perc))), replace = FALSE) 
    sMI <- .normalizeQuantiles2(M[ot=='I', colsel])
    sMII <- .normalizeQuantiles2(M[ot=='II', colsel])
    sUI <- .normalizeQuantiles2(U[ot=='I', colsel])
    sUII <- .normalizeQuantiles2(U[ot=='II', colsel])
    mquan <- list(quantiles = rep(0, nrow(M)),
                    inter = rep(0, nrow(M)),
                    onetwo = ot
                )
    mquan[['quantiles']][ot == 'I'] <- sMI[[1]]
    mquan[['inter']][ot == 'I'] <- sMI[[2]]
    mquan[['quantiles']][ot == 'II'] <- sMII[[1]]
    mquan[['inter']][ot == 'II'] <- sMII[[2]]
    uquan <- list(quantiles = rep(0, nrow(M)),
                    inter = rep(0, nrow(M)),
                    onetwo = ot)
    uquan[['quantiles']][ot == 'I'] <- sUI[[1]]
    uquan[['inter']][ot == 'I'] <- sUI[[2]]
    uquan[['quantiles']][ot == 'II'] <- sUII[[1]]
    uquan[['inter']][ot == 'II'] <- sUII[[2]]
    mquan[['rn']] <- uquan[['rn']] <- rownames(M)

    nmet <- .impose(getMeth(referenceMset), mquan)
    nume <- .impose(getUnmeth(referenceMset), uquan)
    rm(referenceRGset)

    # Everything else continues as normal.
    referenceMset <- minfi::MethylSet(Meth=na.omit(nmet), Unmeth=na.omit(nume), colData=referencePd, annotation(referenceMset))

    if(verbose) message("[estimateCellCounts] Picking probes for composition estimation.\n")
    compData <- minfi:::pickCompProbes(referenceMset, cellTypes = cellTypes, compositeCellType = compositeCellType, probeSelect = probeSelect)
    coefs <- compData$coefEsts
    rm(referenceMset)

    if(verbose) message("[estimateCellCounts] Estimating composition.\n")
    coefdat <- bn[rownames(coefs),]
    rownames(coefdat) <- rownames(coefs)
    counts <- minfi:::projectCellType(coefdat, coefs)
    
    # calculate deconvolution error
    getErrorPerSample = function(applyIndex,
                                 predictedIN = counts,
                                 coefDataIN = coefs,
                                 betasBulkIN = coefdat){
      
      trueBulk = matrix(ncol = 1, nrow = nrow(coefDataIN), data = 0)
      
      RMSE = function(m, o){
        sqrt(mean((m - o)^2))
      }
      
      for (i in 1:ncol(coefDataIN)){
        
        trueBulk[,1] = trueBulk[,1] + coefDataIN[,i]*predictedIN[applyIndex,i]
      }
      
      betasBulkIN = t(apply(betasBulkIN, 1, function(x){x[is.na(x)] = 0; return(x)}))
      
      error = RMSE(trueBulk, betasBulkIN[,applyIndex])
      return(error)
    }
    CellTypePredictionError = sapply(1:nrow(counts), getErrorPerSample)
    counts = cbind(counts, CellTypePredictionError)
    

    if (meanPlot) {
        smeans <- compData$sampleMeans
        smeans <- smeans[order(names(smeans))]
        sampleMeans <- c(colMeans(bn[rownames(coefs),]), smeans)
        sampleColors <- c(rep(1, ncol(coefdat)), 1 + as.numeric(factor(names(smeans))))
        plot(sampleMeans, pch = 21, bg = sampleColors)
        legend("bottomleft", c("blood", levels(factor(names(smeans)))),
               col = 1:7, pch = 15)
    }
    if(returnAll) {
        list(counts = counts, compTable = compData$compTable)
    } else {
        counts
    }
}
