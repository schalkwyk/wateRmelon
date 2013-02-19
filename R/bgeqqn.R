daten1 <-
function(mn, un, onetwo, fudge=100,  ret2=FALSE, ...){

   mnc <- dfsfit(mn, onetwo, ...)
   unc <- dfsfit(un, onetwo, roco=NULL)
   normalizeQuantiles(mnc) -> mncn
   normalizeQuantiles(unc) -> uncn
   beta <- mncn/(mncn + uncn + fudge)
   if (ret2) return (list(methylated=mncn, unmethylated=uncn, beta=beta))
   beta

}
