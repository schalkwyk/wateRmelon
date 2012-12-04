gcoms <-
function(x){ # gets a matrix of ss
   rowSums(x[1:3,], na.rm=TRUE) -> gcoss
   rowSums(x[4:6,], na.rm=TRUE) -> n
   gcoss/n 
}
