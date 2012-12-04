genkme <-
function(y, peaks=c(.2,.5,.8)){ # takes a row of SNP data
   cl <-try( kmeans(as.numeric(y), peaks ), silent=TRUE )
   if (inherits(cl, 'try-error')) {rep(NA,6)}
   else{ c(cl$withinss, cl$size) }
}
