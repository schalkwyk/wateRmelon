nasen <-
function ( mns, uns, onetwo, fudge=100) {

   mns[onetwo=='I' ,] <- normalizeQuantiles(mns[onetwo=='I', ])
   mns[onetwo=='II',] <- normalizeQuantiles(mns[onetwo=='II',])
   
   uns[onetwo=='I' ,] <- normalizeQuantiles(uns[onetwo=='I', ])
   uns[onetwo=='II',] <- normalizeQuantiles(uns[onetwo=='II',])

   mns/(mns + uns + fudge)
}
