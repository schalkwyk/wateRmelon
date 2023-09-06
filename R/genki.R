# LS 27 Jul 2012



#' SNP derived performance metrics for Illumina 450K DNA methylation arrays.
#' 
#' A very simple genotype calling by one-dimensional K-means clustering is
#' performed on each SNP, and for those SNPs where there are three genotypes,
#' the squared deviations are summed for each genotype (similar to a standard
#' deviation for each of allele A homozygote, heterozygote and allele B
#' homozygote).  By default these are further divided by the square root of the
#' number of samples to get a standard error-like statistic.
#' 
#' %% ~~ If necessary, more details than the description above ~~ There are 65
#' well-behaved SNP genotyping probes included on the array.  These each
#' produce a distribution of betas with tight peaks for the three possible
#' genotypes, which will be broadened by technical variation between samples.
#' The spread of the peaks is thus usable as a performance metric.
#' 
#' @param bn a matrix of beta values(default method), a \code{MethyLumiSet}
#' object (\code{methylumi} package), a \code{MethylSet} or \code{RGChannelSet}
#' object (\code{minfi} package) or a \code{exprmethy450} object (\code{IMA}
#' package).
#' @param g vector of SNP names
#' @param se TRUE or FALSE specifies whether to calculate the standard
#' error-like statistic
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... a vector of 3 values for the dispersion of the three
#' genotype peaks (AA, AB, BB : low, medium and high beta values)
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>

#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted)
#' @examples
#' 
#'   #MethyLumiSet method
#'      data(melon)
#'      genki(melon)
#' 
#'   #MethyLumiSet method after normalization
#'      melon.dasen <- dasen(melon)
#'      genki(melon.dasen)
#' 
#' @export genki
genki <- function(bn, g=getsnp(rownames(bn)), se=TRUE ){ # gets a beta matrix
   x     <- apply(bn[g,],1, genkme)
   gcoss <- rowSums(x[1:3,], na.rm=TRUE) 
   n     <- rowSums(x[4:6,], na.rm=TRUE) 
   gcoms <- gcoss/n 
   if( ! se ){return (gcoms) }
   gcoms/sqrt(ncol(bn))
} # returns 3 s{e|m} values
