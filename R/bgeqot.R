dasen <-
function(mns, uns, onetwo, fudge=100, ...){
   mnsc <- dfsfit(mns,  onetwo, ...)  
   unsc <- dfsfit(uns,  onetwo, roco=NULL)
   mnsc[onetwo=='I' ,] <- normalizeQuantiles(mnsc[onetwo=='I', ])
   mnsc[onetwo=='II',] <- normalizeQuantiles(mnsc[onetwo=='II',])

   unsc[onetwo=='I' ,] <- normalizeQuantiles(unsc[onetwo=='I', ])
   unsc[onetwo=='II',] <- normalizeQuantiles(unsc[onetwo=='II',])

   mnsc/( mnsc + unsc + fudge )
}
