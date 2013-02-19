dasen <-
function(mns, uns, onetwo, fudge=100, ret2=FALSE, ...){
   mnsc <- dfsfit(mns,  onetwo, ...)  
   unsc <- dfsfit(uns,  onetwo, roco=NULL)
   mnsc[onetwo=='I' ,] <- normalizeQuantiles(mnsc[onetwo=='I', ])
   mnsc[onetwo=='II',] <- normalizeQuantiles(mnsc[onetwo=='II',])

   unsc[onetwo=='I' ,] <- normalizeQuantiles(unsc[onetwo=='I', ])
   unsc[onetwo=='II',] <- normalizeQuantiles(unsc[onetwo=='II',])
   beta <- mnsc/( mnsc + unsc + fudge )
   if (ret2) return (list(methylated=mnsc,unmethylated=unsc,beta=beta))
   beta
}
