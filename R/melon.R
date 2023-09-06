#' Calculate a full set of 450K normalization/performance metrics
#' 
#' Calculate X-chromosome, SNP and imprinting DMR metrics for a matrix of betas
#' from an Illumina 450K Human DNA methylation array.  Requires precalculated
#' t-test p-values for sex differences, a list of X-chromosome features and of
#' imprinting DMR features.
#' 
#' 
#' @param betas a matrix of betas, each row representing a probe, each column a
#' sample
#' @param pv a vector of p-values such as produced by \code{sextest}, one per
#' row of betas
#' @param X a logical vector of the same length as \code{pv}, indicating
#' whether each probe is mapped to the X-chromosome
#' @param idmr a character vector of probe names known to be in imprinting
#' DMRs.  Can be obtained with \code{iDMR()} or \code{data(iDMR)}
#' @param subset index or character vector giving a subset of betas to be
#' tested
#' @return
#' @returnItem dmrse_row see \code{dmrse_row}
#' @returnItem dmrse_col see \code{dmrse_col}
#' @returnItem dmrse see \code{dmrse }
#' @returnItem gcoms_a see \code{genki }
#' @returnItem gcose_a see \code{genki }
#' @returnItem gcoms_b see \code{genki }
#' @returnItem gcose_b see \code{genki }
#' @returnItem gcoms_c see \code{genki }
#' @returnItem gcose_c see \code{genki }
#' @returnItem seabird see \code{seabi }
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>
#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted)
#' @examples
#' 
#' data(melon)
#' melon.dasen <- dasen(melon)
#' bn <-betas(melon.dasen)
#' X <- melon.dasen@featureData@data$CHR=='X'
#' data(iDMR)
#' sex <- pData(melon.dasen)$sex
#' pv <- sextest(bn,sex)
#' melon.metrics <- metrics(bn, pv, X, idmr = iDMR, subset = NULL) 
#' 
#' @export metrics
metrics <-
function(betas, pv, X, idmr=iDMR, subset=NULL){

   if (! is.null(subset) ) {
      betas<- betas[subset,]
      pv   <- pv   [subset ]
      X    <- X    [subset ]
      
   }
   x        <- genkus(betas)
   gcom     <- gcoms     (x)
   gcos     <- gcose     (x)

   list(

      dmrse_row = dmrse_row (betas, idmr), # 1  formerly SDC - between arrays
      dmrse_col = dmrse_col (betas, idmr), # 2  formerly SDR - between probes
      dmrse     = dmrse     (betas, idmr), # 3  formerly SDO - both
      gcoms_a   = gcom[1]                , # 4  SD like measure low betas 
      gcose_a   = gcos[1]                , # 5  SE like measure low betas
      gcoms_b   = gcom[2]                , # 6  SD like measure med betas
      gcose_b   = gcos[2]                , # 7  SE like measure med betas
      gcoms_c   = gcom[3]                , # 8  SD like measure hi  betas
      gcose_c   = gcos[3]                , # 9  SE like measure hi  betas
      seabird   = 1- seabird(pv, stop=1, X)# 10 1- sexdiff pvalue ROC AUC
  )
}
