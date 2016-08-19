naten <-
function (mn, un, fudge=100, ret2=FALSE ,...) {
   mnn <- normalizeQuantiles(mn)
   unn <- normalizeQuantiles(un)
   beta <- mnn/(mnn + unn + fudge)
   if (ret2) return (list(methylated=mnn,unmethylated=unn,beta=beta))
   beta
}
