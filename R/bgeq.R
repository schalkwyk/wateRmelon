danen <-
function ( mns, uns, onetwo, fudge=100, ...){
   mnsc <- dfsfit(mns,  onetwo, ...)  
   unsc <- dfsfit(uns,  onetwo, roco=NULL)
   mnsc/( mnsc + unsc + fudge )
}
