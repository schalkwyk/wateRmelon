# function to scale away possible sytematic differences between oxybs and bs 450k array data
# that has already been normalized separately eg with dasen

# plan is to 
# 1) calculate rowMeans of betas
# 2) take BS-oxybs
# 3) selct subset with abs(diff) < threshold 
#        threshold is a parameter expressed as a percentile of abs(diff)
# prototype takes 2 methylumi objects and a threshold

oxyscale <- function ( bs, oxybs, threshold=20){

   rmb <- rowMeans(betas(bs))
   omb <- rowMeans(betas(oxybs))
   dif <- rmb - omb 
   # default thresh of 50 gives .25 quantile of negative diff
   thr <- quantile( dif[dif< 0], 1-threshold/200 )
   # equivalent to .25 - .75 quantile of signed diffs thus 50% of data
   use <- abs(dif) < abs(thr)
   cat("using", sum(use), "probes\n")
   # scaling factor:  how much greater is bs than oxybs
   sca <- mean((rmb/omb)[use], na.rm=TRUE )
   cat("scale factor", sca, "\n")
   # multiply oxybs by this
   betas(oxybs) <- betas(oxybs) * sca
   oxybs

}
