sort_order <- function(d, tie=TRUE){
        ## obtain the sorted values and their index
        Si <- sort(d, method = "quick", index.return = TRUE) # NA will be ignored or removed.
        if (tie){
            Si$ix <- NA
        }
        # deal with NA in input d
        nobsj <- length(Si$x)
        #n_1 <- length(d)
        isna <- is.na(d)
        if (sum(isna) > 0) {
            #i <- (0:(n_1 - 1))/(n_1 - 1)
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

uSexQN <- function(mns, uns, ot, chr, cores=1, ret2=FALSE, fudge=100,...){
    stopifnot(nrow(mns) == length(chr))
    stopifnot(nrow(uns) == length(chr))
    stopifnot(nrow(mns) == length(ot))
    stopifnot(nrow(uns) == length(ot))
    stopifnot(length(chr) == length(ot))

    if(!is.logical(chr)){
        is_sex <- grepl('(X|chrX|Y|chrY|23|24)', as.character(chr))
    } else {
        is_sex <- chr
    }

    mns[ot == 'I' , ] <- uSexQNengine(A = mns[ot == 'I' , ], is_sex = is_sex[ot == 'I' ], cores = cores)
    mns[ot == 'II', ] <- uSexQNengine(A = mns[ot == 'II', ], is_sex = is_sex[ot == 'II'], cores = cores)
    uns[ot == 'I' , ] <- uSexQNengine(A = uns[ot == 'I' , ], is_sex = is_sex[ot == 'I' ], cores = cores)
    uns[ot == 'II', ] <- uSexQNengine(A = uns[ot == 'II', ], is_sex = is_sex[ot == 'II'], cores = cores)

    betas = (mns) / (mns+uns+fudge)
    if(ret2){
        return(list(betas = betas, methylated = mns, unmethylated = uns))
    } else {
        return(betas)
    }
}