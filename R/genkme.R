#' Internal functions for genotype-based normalization metrics
#' 
#' genkme - genotype calling with 1d k-means
#' 
#' genkus - apply genkme to available SNPs
#' 
#' getsnp - grep the rs-numbered probes
#' 
#' gcose - calculate between-sample SNP standard error
#' 
#' gcoms - calculate between-sample SNP mean-squared deviation
#' 
#' see \code{\link{genki}}
#' 
#' @aliases genkme genkus getsnp gcose gcoms
#' @param y a vector or matrix of numeric values (betas, between 0 and 1)
#' @param peaks initial values for cluster positions
#' @return
#' 
#' see \code{\link{genki}}
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>
#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted)
#' @export genkme
genkme <-
function(y, peaks=c(.2,.5,.8)){ # takes a row of SNP data
   cl <-try( kmeans(as.numeric(y), peaks ), silent=TRUE )
   if (inherits(cl, 'try-error')) {rep(NA,6)}
   else{ c(cl$withinss, cl$size) }
}
