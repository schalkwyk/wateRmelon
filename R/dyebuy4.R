danet <-
function(mn, un, onetwo, fudge=100, ...){
   mnc <- dfsfit(mn, onetwo, ...)
   unc <- dfsfit(un, onetwo, ...)
   a <- db1(mnc,unc)
   a[[1]]/(a[[1]] + a[[2]] + fudge)
}
