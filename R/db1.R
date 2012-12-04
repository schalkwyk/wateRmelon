db1 <-
function(mn, un ){
   stopifnot(dim(un) == dim(mn))
   a <- dim(un)[2]
   mun <- normalizeQuantiles(cbind(mn,un))
   mnn <- mun[,1:a]
   unn <- mun[,(a+1):(2*a)]
   list(mnn,unn)
}
