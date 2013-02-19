nanes <-
function(mns, uns, onetwo, fudge=100, ret2=FALSE, ...){

   
   mns[onetwo=='I' ,] <- normalizeQuantiles(mns[onetwo=='I', ])
   uns[onetwo=='I' ,] <- normalizeQuantiles(uns[onetwo=='I', ])
   
   a <- db1(mns[onetwo=='II',], uns[onetwo=='II',])
   
   mns[onetwo=='II' ,] <- a[[1]]
   uns[onetwo=='II' ,] <- a[[2]]
   
   beta <- mns/(mns + uns + fudge)
   if (ret2) return (list(methylated=mns,unmethylated=uns,beta=beta))
   beta

}
