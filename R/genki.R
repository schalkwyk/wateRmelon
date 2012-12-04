# LS 27 Jul 2012

genki <- function(bn, g=getsnp(rownames(bn)), se=TRUE ){ # gets a beta matrix
   x     <- apply(bn[g,],1, genkme)
   gcoss <- rowSums(x[1:3,], na.rm=TRUE) 
   n     <- rowSums(x[4:6,], na.rm=TRUE) 
   gcoms <- gcoss/n 
   if( ! se ){return (gcoms) }
   gcoms/sqrt(ncol(bn))
} # returns 3 s{e|m} values
