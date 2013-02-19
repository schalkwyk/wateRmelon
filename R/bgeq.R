danen <-
function ( mns, uns, onetwo, fudge=100, ret2=FALSE, ...){
   mnsc <- dfsfit(mns,  onetwo, ...)  
   unsc <- dfsfit(uns,  onetwo, roco=NULL)
   beta <- mnsc/( mnsc + unsc + fudge )
   if (ret2) return (list(methylated=mnsc,unmethylated=unsc,beta=beta))
   beta
}
