nanes <-
function(mns, uns, onetwo, fudge=100){

   
   mns[onetwo=='I' ,] <- normalizeQuantiles(mns[onetwo=='I', ])
   uns[onetwo=='I' ,] <- normalizeQuantiles(uns[onetwo=='I', ])
   
   a <- db1(mns[onetwo=='II',], uns[onetwo=='II',])
   
   mns[onetwo=='II' ,] <- a[[1]]
   uns[onetwo=='II' ,] <- a[[2]]
   
   mns/(mns + uns + fudge)

}
