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
