nasen <-
function ( mns, uns, onetwo, ret2=FALSE, fudge=100,...) {

   mns[onetwo=='I' ,] <- normalizeQuantiles(mns[onetwo=='I', ])
   mns[onetwo=='II',] <- normalizeQuantiles(mns[onetwo=='II',])
   
   uns[onetwo=='I' ,] <- normalizeQuantiles(uns[onetwo=='I', ])
   uns[onetwo=='II',] <- normalizeQuantiles(uns[onetwo=='II',])

   beta <-  mns/(mns + uns + fudge)
   if (ret2) return (list(methylated=mns, unmethylated=uns, beta=beta))
   beta
}
