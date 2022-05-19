#' 
#' adjustedDasen
#' 
#' @description adjustedDasen utilizes dasen normliasation to normalise autosomal
#' CpGs, and infers the sex chromosome linked CpGs by linear interpolation on 
#' corrected autosomal CpGs.
#'
#' @param mns matrix of methylated signal intensities, samples in column and 
#' probes in row.
#' @param uns matrix of unmethylated signal intensities, samples in column and 
#' probes in row.
#' @param onetwo character vector or factor of length nrow(mns) indicating assay
#'  type 'I' or 'II'.
#' @param chr character vector stores the mapped chromosomes for all probes, e.g. 
#' chr <- c('1', 'X', '21', ..., 'Y').
#' @param offset_fit logical (default is TRUE). To use dasen, set it TRUE; to use
#' nasen, set it FALSE.
#' @param cores an integer(e.g. 8) defines the number of cores to parallel processing. 
#' Default value is 1, set to -1 to use all available cores.
#' @param ret2 logical (default is FALSE), if TRUE, returns a list of intensities 
#' and betas instead of a naked matrix of betas.
#' @param fudge default 100, a value added to total intensity to prevent denominators 
#' close to zero when calculating betas, e.g. betas <- mns / (mns + uns + fudge).
#' @param ... additional argument roco for dfsfit giving Sentrix rows and 
#' columns.  This allows a background gradient model to be fit. This is split 
#' from data column names by default.  roco=NULL disables model fitting (and 
#' speeds up processing), otherwise roco can be supplied as a character vector 
#' of strings like 'R01C01' (only 3rd and 6th characters used).
#'
#' @return a matrix of normalised beta values.
#' @export
#' adjustedDasen
#'
#' @references 
#' A data-driven approach to preprocessing Illumina 450K methylation array data, 
#' Pidsley et al, BMC Genomics. \cr
#' interpolatedXY: a two-step strategy to normalise DNA methylation 
#' microarray data avoiding sex bias, Wang et al., 2021.
#' @examples
#' data(melon)
#' normalised_betas <- adjustedDasen(mns = methylated(melon), uns = unmethylated(melon), onetwo = fData(melon)[,fot(melon)], chr = fData(melon)$CHR, cores=1)
#' ## if input is an object of methylumiset or methylset
#' normalised_betas <- adjustedDasen(melon)
#' 
adjustedDasen <- function(mns, uns, onetwo, chr, offset_fit=TRUE, cores=1, ret2=FALSE, fudge=100,...){
    stopifnot(nrow(mns) == length(chr))
    stopifnot(nrow(uns) == length(chr))
    stopifnot(nrow(mns) == length(onetwo))
    stopifnot(nrow(uns) == length(onetwo))
    stopifnot(length(chr) == length(onetwo))
    
    if(!is.logical(chr)){
        is_sex <- grepl('(X|chrX|Y|chrY|23|24)', as.character(chr))
    } else {
        is_sex <- chr
    }    
    
    if (cores < 1) {
        cores <- detectCores()
    }
    
    if(Sys.info()["sysname"] != "Linux"){
        cores <- 1
    }
    
    ## to use 'nasen', set offset_fit=FALSE
    if(offset_fit){        
        mns <- p_dfsfit(mns,  onetwo, cores=cores)
        uns <- p_dfsfit(uns,  onetwo, roco=NULL, cores=cores)    
    }
    
    mns[onetwo == 'I' , ] <- uSexQNengine(A = mns[onetwo == 'I' , ], is_sex = is_sex[onetwo == 'I' ], cores = cores)
    mns[onetwo == 'II', ] <- uSexQNengine(A = mns[onetwo == 'II', ], is_sex = is_sex[onetwo == 'II'], cores = cores)
    uns[onetwo == 'I' , ] <- uSexQNengine(A = uns[onetwo == 'I' , ], is_sex = is_sex[onetwo == 'I' ], cores = cores)
    uns[onetwo == 'II', ] <- uSexQNengine(A = uns[onetwo == 'II', ], is_sex = is_sex[onetwo == 'II'], cores = cores)
    
    betas = (mns) / (mns+uns+fudge)
    if(ret2){
        return(list(betas = betas, methylated = mns, unmethylated = uns))
    } else {
        return(betas)
    }
}


sort_order <- function(d, tie=TRUE){
    ## obtain the sorted values and their index
    Si <- sort(d, method = "quick", index.return = TRUE) # NA will be ignored or removed.
    if (tie){
        Si$ix <- NA
    }
    # deal with NA in input d
    nobsj <- length(Si$x)
    n_1 <- length(d)
    isna <- is.na(d)
    if (sum(isna) > 0) {
        i <- (0:(n_1 - 1))/(n_1 - 1)
        Si$x <- approx((0:(nobsj - 1))/(nobsj - 1), Si$x, i, ties = list("ordered", mean))$y  # Si$x will not contain NAs any more.
        if (!tie) {
            O_i <- rep(NA, n_1)
            O_i[!isna] <- ((1:n_1)[!isna])[Si$ix]
            Si$ix <- O_i
        }
    }
    return(Si)
}


tie_norm <- function(d, is_sex, rank2mean){
    ## normalise d differently on autosomes and XY
    d_sex <- d[is_sex]
    d_autosome <- d[!is_sex]
    r_autosome <- rank(d_autosome) # NA will be counted and placed at the end.
    # Get the ranks of sexual cpgs based on ranks of  autosomal cpgs;
    # rule=2 means the value at the closest data extreme is used when new x is greater than max(x)
    r_sex <- approx(d_autosome, r_autosome, d_sex, ties = mean, rule=2)$y
    
    # Produce the final values of non-NA  autosomal cpgs based on their ranks
    notna <- !is.na(d_autosome)
    nobsj <- sum(notna)
    d[!is_sex][notna] <- rank2mean((r_autosome[notna] - 1)/(nobsj - 1))
    # Produce the final values of non-NA sexual cpgs based on their ranks
    notna_sex <- !is.na(d_sex)
    d[is_sex][notna_sex] <- rank2mean((r_sex[notna_sex] - 1)/(nobsj - 1))
    return(d)
}


uSexQNengine <- function(A, is_sex, cores=1) {
    ## A: a dataframe or matrix;
    ## chr: a vector, like c('1', '2', 'X', 'Y')
    stopifnot(nrow(A) == length(is_sex))
    A <- data.frame(A, check.names=FALSE)
    A_autosome <- A[!is_sex, ]
    n_1 <- nrow(A_autosome)
    sort_Aa <- mclapply(A_autosome, sort_order, mc.cores=cores, tie=TRUE)
    S_autosome <- sapply(sort_Aa, function(x) x$x)
    m_autosome <- rowMeans(S_autosome)
    
    # Get a function which gives relationships between orders and mean values.
    i <- (0:(n_1 - 1))/(n_1 - 1)
    rank2mean <- approxfun(i, m_autosome, ties = list("ordered", mean))
    
    #rm(S_autosome, A_autosome, sort_Aa)
    # For each sample, find its normalised values
    A <- mclapply(A, tie_norm, is_sex=is_sex, mc.cores=cores, rank2mean=rank2mean)
    A <- sapply(A, function(x) x)
    
    return(A)
}


p_dfsfit <- function (mn, onetwo, cores=1, roco=substring(colnames(mn), regexpr("R0[1-9]C0[1-9]", colnames(mn))), ...){
    mn <- data.frame(mn, check.names=FALSE)
    mdf <- mclapply(mn, dfs2, onetwo, mc.cores=cores)
    mdf <- sapply(mdf, function(x) x)
    if (!is.null(roco)) {
        scol <- as.numeric(substr(roco, 6, 6))
        srow <- as.numeric(substr(roco, 3, 3))
        fit <- try(lm(mdf ~ srow + scol), silent = TRUE)
        if (!inherits(fit, "try-error")) {
            mdf <- fit$fitted.values
        }
        else {
            message("Sentrix position model failed, skipping")
        }
    }
    otcor <- matrix(rep(mdf, sum(onetwo == "I")), byrow = T, nrow = sum(onetwo == "I"))
    mn[onetwo == "I", ] <- mn[onetwo == "I", ] - otcor
    mn
}
