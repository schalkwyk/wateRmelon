#' Calculate ROC area-under-curve for X-chromosome sex differences (internal
#' function for calculating the seabi metric)
#' 
#' This is a wrapper for the prediction and performance functions from the ROCR
#' package that takes a vector of p-values and a vector of true or false for
#' being on the X.  See \code{seabi} function which does everything.
#' 
#' 
#' @param pr a vector of p-values, such as calculated by \code{seabird}
#' @param stop fraction for partial area under curve.  For example 0.1 gives
#' you the area for the lowest 10\% of p-values.
#' @param X logical vector the same length as pv, true for features mapped to
#' X-chromosome
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... Returns an area value between 0 and 1, where 1 is the best
#' possible performance.
#' @author Leonard C Schalkwyk 2012 Leonard.Schalkwyk@@kcl.ac.uk
#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted) %% ~put references to the literature/web site here ~
#' @export seabird
seabird <-
function(pr, stop=1, X){ # sexdiff pvalue ROC AUC
                                    # pr : pvals from sextest
                                    # stop: fraction for partial AUC
                                    # X: logical vector-probe on X?
# ROCR version
#    pr <- prediction(1 - pr, X)
#    unlist(performance(pr, "auc", fpr.stop = stop)@y.values)
# pROC version
#    as.numeric(auc(formula=X~pr, data=NULL)) 
# ROC version (bioconductor)

   pAUC(rocdemo.sca(truth=X, data=1-pr), stop) 

}
   





#' Calculate a performance metric based on male-female differences for Illumina
#' methylation 450K arrays
#' 
#' Calculates an area under ROC curve - based metric for Illumina 450K data
#' using a t-test for male-female difference as the predictor for X-chromosome
#' location of probes.  The metric is 1-area so that small values indicate good
#' performance, to match our other, standard error based metrics
#' \code{\link{gcose}} and \code{\link{dmrse}}. Note that this requires both
#' male and female samples of known sex and can be slow to compute due to
#' running a t-test on every probe.
#' 
#' 
#' @param bn a matrix of betas (default method) or an object containing betas
#' i.e. a \code{MethyLumiSet} object (\code{methylumi} package), a
#' \code{MethylSet} or \code{RGChannelSet} object (\code{minfi} package) or a
#' \code{exprmethy450} object (\code{IMA} package).
#' @param stop partial area under curve is calculated if stop value <1 is
#' provided
#' @param sex a factor giving the sex of each sample (column)
#' @param X a logical vector of length equal to the number of probes, true for
#' features mapped to X-chromosome
#' @return
#' 
#' a value between 0 and 1.  values close to zero indicate high data quality as
#' judged by the ability to discriminate male from female X-chromosome DNA
#' methylation.
#' 
#' %% ~Describe the value returned %% If it is a LIST, use %% \item{comp1
#' }{Description of 'comp1'} %% \item{comp2 }{Description of 'comp2'} %% ...
#' @author leonard.schalkwyk@@kcl.ac.uk
#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted)
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#'    library(methylumi)
#'    data(melon)
#'    sex  <- pData(melon)$sex
#'    X    <- melon@featureData@data$CHR=='X'
#'    seabi(betas(melon), sex=sex, X=X)
#' 
#' # methylumi method
#'    seabi(melon, sex=sex, X=X)
#' 
#' @export seabi
seabi <- 
function (bn, stop=1, sex, X){
   1 - seabird( sextest(bn, sex), stop, X )
}
