nanet <-
function(mn, un, fudge=100, ret2=FALSE,...){
   a <- db1(mn,un)
   beta <- a[[1]]/(a[[1]] + a[[2]] + fudge)
   if (ret2) return (list(methylated=a[[1]],unmethylated=a[[2]],beta=beta))
   beta
}
