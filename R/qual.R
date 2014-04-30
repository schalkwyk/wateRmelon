qual <- function(x,y){ # normalized and original betas
  dif  <- x - y
  rmsd <- sqrt(colMeans(dif^2, na.rm=TRUE))
  sdd  <- apply(dif, 2, sd, na.rm=TRUE) 
  srms <- rmsd/sdd
  data.frame(rmsd,sdd,srms) 
}
