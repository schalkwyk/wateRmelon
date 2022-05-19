
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