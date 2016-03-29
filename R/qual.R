qual <- function(norm,raw){ # normalized and original betas
  dif  <- norm - raw
  rmsd <- sqrt(colMeans(dif^2, na.rm=TRUE))
  sdd  <- apply(dif, 2, sd, na.rm=TRUE) 
  sadd  <- apply(abs(dif), 2, sd, na.rm=TRUE) 
  srms <- rmsd/sdd
  data.frame(rmsd,sdd,sadd,srms) 
}
