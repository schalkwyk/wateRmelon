danet <-
function(mn, un, onetwo, fudge=100, ret2=FALSE, ...){
   mnc <- dfsfit(mn, onetwo, ...)
   unc <- dfsfit(un, onetwo, ...)
   a <- db1(mnc,unc)
   beta <- a[[1]]/(a[[1]] + a[[2]] + fudge)
   if (ret2) return (list(methylated=a[[1]],unmethylated=a[[2]],beta=beta))
   beta
}
