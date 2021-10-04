### ORIGINAL AUTHOR: Jean-Philippe Fortin, Sept 24 2013 (Functional normalization of 450k methylation array data improves replication in large cancer studies, Genome Biology, 2014)
### Adopted by Yucheng Wang

######################################################
## Functional normalization of the 450k array
## Jean-Philippe Fortin
## Sept 24 2013
#####################################################

##


#' adjustedFunnorm
#' 
#' @description adjustedFunnorm utilizes functional normliasation to normalise autosomal
#' CpGs, and infers the sex chromosome linked CpGs by linear interpolation on 
#' corrected autosomal CpGs.
#'
#' @param rgSet An object of class "RGChannelSet".
#' @param nPCs Number of principal components from the control probes PCA.
#' @param sex An optional numeric vector containing the sex of the samples.
#' @param bgCorr Should the NOOB background correction be done, prior to 
#' functional normalization (see "preprocessNoob")
#' @param dyeCorr Should dye normalization be done as part of the NOOB 
#' background correction (see "preprocessNoob")?
#' @param keepCN Should copy number estimates be kept around? Setting to 'FALSE'
#'  will decrease the size of the output object significantly.
#' @param ratioConvert Should we run "ratioConvert", ie. should the output be a 
#' "GenomicRatioSet" or should it be kept as a "GenomicMethylSet"; the latter 
#' is for experts.
#' @param verbose Should the function be verbose?
#' 
#' @return an object of class "GenomicRatioSet", unless "ratioConvert=FALSE" in 
#' which case an object of class "GenomicMethylSet".
#' @export 
#' adjustedFunnorm
#' @references  
#' Functional normalization of 450k methylation array data improves replication 
#' in large cancer studies, Fortin et al., 2014, Genome biology. \cr
#' interpolatedXY: a two-step strategy to normalise DNA methylation 
#' microarray data avoiding sex bias, Wang et al., 2021.
#'
#' @examples
#' \dontrun{
#' GRset <- adjustedFunnorm(RGSet)
#' }
#' 
adjustedFunnorm <- function(rgSet, nPCs=2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE, verbose = TRUE) {
    
    .isMatrixBackedOrStop(rgSet, "adjustedFunnorm")
    
    .isRGOrStop(rgSet)
    rgSet <- updateObject(rgSet) ## FIXM: might not KDH: technically, this should not be needed, but might be nice
    
    # Background correction and dye bias normalization:
    if (bgCorr){
        if(verbose && dyeCorr) {
            message("[adjustedFunnorm] Background and dye bias correction with noob")
        } else {
            message("[adjustedFunnorm] Background correction with noob")
        }
        gmSet <- preprocessNoob(rgSet, dyeCorr = dyeCorr)
        if(verbose) message("[adjustedFunnorm] Mapping to genome")
        gmSet <- mapToGenome(gmSet)
    } else {
        if(verbose) message("[adjustedFunnorm] Mapping to genome")
        gmSet <- mapToGenome(rgSet)
    }
    
    subverbose <- max(as.integer(verbose) - 1L, 0)
    
    if(verbose) message("[adjustedFunnorm] Quantile extraction")
    extractedData <- .extractFromRGSet450k(rgSet)
    rm(rgSet)
    
    if (is.null(sex)) {
        gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
        sex <- rep(1L, length(gmSet$predictedSex))
        sex[gmSet$predictedSex == "F"] <- 2L
    }
    if(verbose) message("[adjustedFunnorm] Normalization")
    if(keepCN) {
        CN <- getCN(gmSet)
    }
    gmSet <- .adjusted_normalizeFunnorm450k(object = gmSet, extractedData = extractedData,
                                            sex = sex, nPCs = nPCs, verbose = subverbose)
    preprocessMethod <- c(preprocessMethod(gmSet),
                          mu.norm = sprintf("Funnorm, nPCs=%s", nPCs))
    if(ratioConvert) {
        grSet <- ratioConvert(gmSet, type = "Illumina", keepCN = keepCN)
        if(keepCN) {
            assay(grSet, "CN") <- CN
        }
        grSet@preprocessMethod <- preprocessMethod
        return(grSet)
    } else {
        gmSet@preprocessMethod <- preprocessMethod
        return(gmSet)
    }
}


.getFunnormIndices <- function(object) {
    ## WYC
    .isGenomicOrStop(object)
    probeType <- getProbeType(object, withColor = TRUE)
    # autosomal <- (seqnames(object) %in% paste0("chr", 1:22))
    indices <- list(IGrn = (probeType == "IGrn"),
                    IRed = (probeType == "IRed"),
                    II = (probeType == "II" ),
                    XY = as.vector(seqnames(object)) %in% c("chrX", "chrY"))
    indices
}


.adjusted_normalizeFunnorm450k <- function(object, extractedData, nPCs, sex, verbose = TRUE) {
    #normalizeQuantiles <- function(matrix, indices, sex = NULL) {
    #    matrix <- matrix[indices,,drop=FALSE]
    #    ## uses probs, model.matrix, nPCS, through scoping)
    #    oldQuantiles <- t(colQuantiles(matrix, probs = probs))
    #    if(is.null(sex)) {
    #        newQuantiles <- .returnFit(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs)
    #    } else {
    #        newQuantiles <- .returnFitBySex(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs, sex = sex)
    #    }
    #    .normalizeMatrix(matrix, newQuantiles)
    #}
    
    interpolatedXY <- function(ra_signal, na_signal, ru_signal){
        # construct a function which reflects relationships between orders and final norm values.
        n_1 <- length(ra_signal)
        rank2mean <- approxfun((0:(n_1 - 1))/(n_1 - 1), sort(na_signal, method = "quick"), ties = list("ordered", mean)) 
        
        rank_autosome <- rank(ra_signal) # NA will be counted and placed at the end.
        # Get the ranks of sexual cpgs based on ranks of  autosomal cpgs;
        # rule=2 means the value at the closest data extreme is used when new x is greater than max(x)
        rank_sex <- approx(ra_signal, rank_autosome, ru_signal, ties = mean, rule=2)$y
        
        # Produce the final values of non-NA sexual cpgs based on their ranks
        notNA <- !is.na(ru_signal)
        ru_signal[notNA] <- rank2mean((rank_sex[notNA] - 1)/(n_1 - 1))
        ru_signal
    } 
    
    unbiased_normalizeQuantiles <- function(mat, indices, sex_probe = NULL) {
        mat_auto <- mat[(indices & !sex_probe),,drop=FALSE]
        mat_sex <- mat[(indices & sex_probe),,drop=FALSE]
        ## uses probs, model.matrix, nPCS, through scoping)
        oldQuantiles <- t(colQuantiles(mat_auto, probs = probs))
        newQuantiles <- .returnFit(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs)
        n_matrix <- .normalizeMatrix(mat_auto, newQuantiles)
        for(j in 1:ncol(n_matrix)){
            mat_sex[, j] <- interpolatedXY(mat_auto[, j], n_matrix[, j], mat_sex[, j])
        }
        mat[(indices & sex_probe), ] <- mat_sex
        mat[(indices & !sex_probe), ] <- n_matrix
        mat
    }
    
    indicesList <- .getFunnormIndices(object)
    model.matrix <- .buildControlMatrix450k(extractedData)
    probs <- seq(from = 0, to = 1, length.out = 500)
    Meth <- getMeth(object)
    Unmeth <- getUnmeth(object)
    if (nPCs > 0){
        for (type in c("IGrn", "IRed", "II")) {
            indices <- indicesList[[type]]
            if(length(indices) > 0) {
                if(verbose) message(sprintf("[InterpolatedXY adjustedFunnorm] Normalization of the %s probes", type))
                Unmeth <- unbiased_normalizeQuantiles(Unmeth, indices=indices, sex_probe=indicesList$XY)
                Meth <- unbiased_normalizeQuantiles(Meth, indices=indices, sex_probe=indicesList$XY)
            }
        }
    }
    assay(object, "Meth") <- Meth
    assay(object, "Unmeth") <- Unmeth
    return(object)
}


### To extract quantiles and control probes from rgSet
.extractFromRGSet450k <- function(rgSet) {
    rgSet <- updateObject(rgSet)
    controlType <- c("BISULFITE CONVERSION I",
                     "BISULFITE CONVERSION II",
                     "EXTENSION",
                     "HYBRIDIZATION",
                     "NEGATIVE",
                     "NON-POLYMORPHIC",
                     "NORM_A",
                     "NORM_C",
                     "NORM_G",
                     "NORM_T",
                     "SPECIFICITY I",
                     "SPECIFICITY II",
                     "TARGET REMOVAL",
                     "STAINING")
    
    array <- annotation(rgSet)[["array"]]
    ## controlAddr <- getControlAddress(rgSet, controlType = controlType, asList = TRUE)
    ctrls <- getProbeInfo(rgSet, type = "Control")
    ctrls <- data.frame(ctrls)
    if(!all(controlType %in% ctrls$Type))
        stop("The `rgSet` does not contain all necessary control probes")
    
    ctrlsList <- split(ctrls, ctrls$Type)[controlType]
    redControls <- getRed(rgSet)[ctrls$Address,,drop=FALSE]
    redControls <- lapply(ctrlsList, function(ctl) redControls[ctl$Address,,drop=FALSE])
    greenControls <- getGreen(rgSet)[ctrls$Address,,drop=FALSE]
    greenControls <- lapply(ctrlsList, function(ctl) greenControls[ctl$Address,,drop=FALSE])
    
    ## Extraction of the undefined negative control probes
    oobRaw <- getOOB(rgSet)
    probs <- c(0.01, 0.50, 0.99)
    greenOOB <- t(colQuantiles(oobRaw$Grn, na.rm = TRUE, probs = probs))
    redOOB   <- t(colQuantiles(oobRaw$Red, na.rm=TRUE,  probs = probs))
    oob      <- list(greenOOB = greenOOB, redOOB = redOOB)
    
    return(list(
        greenControls = greenControls,
        redControls = redControls,
        oob = oob, ctrlsList = ctrlsList,
        array = array))
}


## Extraction of the Control matrix
.buildControlMatrix450k <- function(extractedData) {
    getCtrlsAddr <- function(exType, index) {
        ctrls <- ctrlsList[[index]]
        addr <- ctrls$Address
        names(addr) <- ctrls$ExtendedType
        na.omit(addr[exType])
    }
    
    array <- extractedData$array
    greenControls <- extractedData$greenControls
    redControls <- extractedData$redControls
    controlNames <- names(greenControls)
    ctrlsList <- extractedData$ctrlsList
    
    ## Bisulfite conversion extraction for probe type II:
    index <- match("BISULFITE CONVERSION II", controlNames)
    redControls.current <- redControls[[ index ]]
    bisulfite2 <- colMeans2(redControls.current, na.rm = TRUE)
    
    ## Bisulfite conversion extraction for probe type I:
    index <- match("BISULFITE CONVERSION I", controlNames)
    if (array=="IlluminaHumanMethylation450k"){
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3), index = index)
    } else {
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", c("-", "-"), 1:2), index = index)
    }
    greenControls.current <- greenControls[[ index ]][addr,,drop=FALSE]
    if (array=="IlluminaHumanMethylation450k"){
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 4:6), index = index)
    } else {
        addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 3:5), index = index)
    }
    redControls.current <- redControls[[ index ]][addr,, drop=FALSE]
    if (nrow(redControls.current)==nrow(greenControls.current)){
        bisulfite1 <- colMeans2(redControls.current + greenControls.current, na.rm = TRUE)
    } else {
        bisulfite1 <- colMeans2(redControls.current, na.rm=TRUE) + colMeans2(greenControls.current, na.rm = TRUE)
    }
    
    
    ## Staining
    index <- match("STAINING", controlNames)
    addr <- getCtrlsAddr(exType = "Biotin (High)", index = index)
    stain.green <- t(greenControls[[ index ]][addr,,drop=FALSE])
    addr <- getCtrlsAddr(exType = "DNP (High)", index = index)
    stain.red <- t(redControls[[ index ]][addr,, drop=FALSE ])
    
    ## Extension
    index <-    match("EXTENSION", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("A", "T")), index = index)
    extension.red <- t(redControls[[index]][addr,,drop=FALSE])
    colnames(extension.red) <- paste0("extRed", 1:ncol(extension.red))
    addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("C", "G")), index = index)
    extension.green <- t(greenControls[[index]][addr,,drop=FALSE])
    colnames(extension.green) <- paste0("extGrn", 1:ncol(extension.green))
    
    ## Hybridization should be monitored only in the green channel
    index <- match("HYBRIDIZATION", controlNames)
    hybe <- t(greenControls[[index]])
    colnames(hybe) <- paste0("hybe", 1:ncol(hybe))
    
    ## Target removal should be low compared to hybridization probes
    index <- match("TARGET REMOVAL", controlNames)
    targetrem <- t(greenControls[[index]])
    colnames(targetrem) <- paste0("targetrem", 1:ncol(targetrem))
    
    ## Non-polymorphic probes
    index <- match("NON-POLYMORPHIC", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("A", "T")), index = index)
    nonpoly.red <- t(redControls[[index]][addr, ,drop=FALSE])
    colnames(nonpoly.red) <- paste0("nonpolyRed", 1:ncol(nonpoly.red))
    addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("C", "G")), index = index)
    nonpoly.green <- t(greenControls[[index]][addr, ,drop=FALSE])
    colnames(nonpoly.green) <- paste0("nonpolyGrn", 1:ncol(nonpoly.green))
    
    ## Specificity II
    index <- match("SPECIFICITY II", controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    spec2.green <- t(greenControls.current)
    colnames(spec2.green) <- paste0("spec2Grn", 1:ncol(spec2.green))
    spec2.red <- t(redControls.current)
    colnames(spec2.red) <- paste0("spec2Red", 1:ncol(spec2.red))
    spec2.ratio <- colMeans2(greenControls.current, na.rm = TRUE) /
        colMeans2(redControls.current, na.rm = TRUE)
    
    ## Specificity I
    index <- match("SPECIFICITY I", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 1:3), index = index)
    greenControls.current <- greenControls[[index]][addr,,drop=FALSE]
    redControls.current <- redControls[[index]][addr,,drop=FALSE]
    spec1.green <- t(greenControls.current)
    colnames(spec1.green) <- paste0("spec1Grn", 1:ncol(spec1.green))
    spec1.ratio1 <- colMeans2(redControls.current, na.rm = TRUE) /
        colMeans2(greenControls.current, na.rm = TRUE)
    
    index <- match("SPECIFICITY I", controlNames) # Added that line
    addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 4:6), index = index)
    greenControls.current <- greenControls[[index]][addr,,drop=FALSE]
    redControls.current <- redControls[[index]][addr,,drop=FALSE]
    spec1.red <- t(redControls.current)
    colnames(spec1.red) <- paste0("spec1Red", 1:ncol(spec1.red))
    spec1.ratio2 <- colMeans2(greenControls.current, na.rm = TRUE) /
        colMeans2(redControls.current, na.rm = TRUE)
    spec1.ratio <- (spec1.ratio1 + spec1.ratio2) / 2
    
    ## Normalization probes:
    index <- match(c("NORM_A"), controlNames)
    normA <- colMeans2(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_T"), controlNames)
    normT <- colMeans2(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_C"), controlNames)
    normC <- colMeans2(greenControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_G"), controlNames)
    normG <- colMeans2(greenControls[[index]], na.rm = TRUE)
    
    dyebias <- (normC + normG)/(normA + normT)
    
    oobG <- extractedData$oob$greenOOB
    oobR <- extractedData$oob$redOOB
    oob.ratio <- oobG[2,]/oobR[2,]
    oobG <- t(oobG)
    colnames(oobG) <- paste0("oob", c(1,50,99))
    
    model.matrix <- cbind(
        bisulfite1, bisulfite2, extension.green, extension.red, hybe,
        stain.green, stain.red, nonpoly.green, nonpoly.red,
        targetrem, spec1.green, spec1.red, spec2.green, spec2.red, spec1.ratio1,
        spec1.ratio, spec2.ratio, spec1.ratio2, normA, normC, normT, normG, dyebias,
        oobG, oob.ratio)
    
    
    ## Imputation
    for (colindex in 1:ncol(model.matrix)) {
        if(any(is.na(model.matrix[,colindex]))) {
            column <- model.matrix[,colindex]
            column[is.na(column)]    <- mean(column, na.rm = TRUE)
            model.matrix[ , colindex] <- column
        }
    }
    
    ## Scaling
    model.matrix <- scale(model.matrix)
    
    ## Fixing outliers
    model.matrix[model.matrix > 3] <- 3
    model.matrix[model.matrix < (-3)] <- -3
    
    ## Rescaling
    model.matrix <- scale(model.matrix)
    
    return(model.matrix)
}


### Return the normalized quantile functions
.returnFit <- function(controlMatrix, quantiles, nPCs) {
    stopifnot(is.matrix(quantiles))
    stopifnot(is.matrix(controlMatrix))
    stopifnot(ncol(quantiles) == nrow(controlMatrix))
    ## Fixing potential problems with extreme quantiles
    quantiles[1,] <- 0
    quantiles[nrow(quantiles),] <- quantiles[nrow(quantiles) - 1,] + 1000
    meanFunction <- rowMeans2(quantiles)
    res <- quantiles - meanFunction
    controlPCs <- prcomp(controlMatrix)$x[,1:nPCs,drop=FALSE]
    design <- model.matrix(~controlPCs)
    fits <- lm.fit(x = design, y = t(res))
    newQuantiles <- meanFunction + t(fits$residuals)
    newQuantiles <- .regularizeQuantiles(newQuantiles)
    return(newQuantiles)
}

.returnFitBySex <- function(controlMatrix, quantiles, nPCs, sex) {
    stopifnot(is.matrix(quantiles))
    stopifnot(is.matrix(controlMatrix))
    stopifnot(ncol(quantiles) == nrow(controlMatrix))
    sex    <- as.character(sex)
    levels <- unique(sex)
    nSexes <- length(levels)
    if (nSexes == 2) {
        sex1 <- sum(sex == levels[1])
        sex2 <- sum(sex == levels[2])
        
    } else {
        sex1 <- sum(sex == levels[1])
        sex2 <- 0
    }
    
    ## When normalization should not be performed by sex separately:
    if ((sex1 <= 10) | (sex2 <= 10)) {
        newQuantiles <- .returnFit(controlMatrix = controlMatrix,
                                   quantiles = quantiles,
                                   nPCs = nPCs)
    } else {
        quantiles1 <- quantiles[, sex == levels[1]]
        controlMatrix1 <- controlMatrix[sex == levels[1], ]
        
        newQuantiles1 <- .returnFit(controlMatrix = controlMatrix1,
                                    quantiles = quantiles1,
                                    nPCs = nPCs)
        
        quantiles2 <- quantiles[, sex == levels[2]]
        controlMatrix2 <- controlMatrix[sex == levels[2], ]
        
        newQuantiles2 <- .returnFit(controlMatrix = controlMatrix2,
                                    quantiles = quantiles2,
                                    nPCs = nPCs)
        
        newQuantiles <- quantiles
        newQuantiles[, sex == levels[1]] <- newQuantiles1
        newQuantiles[, sex == levels[2]] <- newQuantiles2
    }
    
    return(newQuantiles)
}


### Normalize a matrix of intensities
.normalizeMatrix <- function(intMatrix, newQuantiles) {
    ## normMatrix <- matrix(NA, nrow(intMatrix), ncol(intMatrix))
    n <- nrow(newQuantiles)
    normMatrix <- sapply(1:ncol(intMatrix), function(i) {
        crtColumn <- intMatrix[ , i]
        crtColumn.reduced <- crtColumn[!is.na(crtColumn)]
        ## Generation of the corrected intensities:
        target <- sapply(1:(n-1), function(j) {
            start <- newQuantiles[j,i]
            end <- newQuantiles[j+1,i]
            if (!isTRUE(all.equal(start,end))){
                sequence <- seq(start, end,( end-start)/n)[-(n+1)]
            } else {
                sequence <- rep(start, n)
            }
            return(sequence)
        })
        target <- as.vector(target)
        result <- preprocessCore::normalize.quantiles.use.target(matrix(crtColumn.reduced), target)
        return(result)
    })
    return(normMatrix)
}

# To ensure a monotonically increasing and non-negative quantile function
# Necessary for pathological cases
.regularizeQuantiles <- function(x){
    x[x<0] <- 0
    colCummaxs(x)
}


## WYC
.isMatrixBackedOrStop <- function(object, FUN) {
    if (!.isMatrixBacked(object)) {
        stop("'", FUN, "()' only supports matrix-backed minfi objects.",
             call. = FALSE)
    }
}


.isMatrixBacked <- function(object) {
    stopifnot(is(object, "SummarizedExperiment"))
    all(vapply(assays(object), is.matrix, logical(1L)))
}


.isRGOrStop <- function(object) {
    if (!is(object, "RGChannelSet")) {
        stop("object is of class '", class(object), "', but needs to be of ",
             "class 'RGChannelSet' or 'RGChannelSetExtended'")
    }
}


.isGenomicOrStop <- function(object) {
    if (!is(object, "GenomicMethylSet") && !is(object, "GenomicRatioSet")) {
        stop("object is of class '", class(object), "', but needs to be of ",
             "class 'GenomicMethylSet' or 'GenomicRatioSet'")
    }
}

