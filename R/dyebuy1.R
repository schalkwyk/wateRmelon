nanet <-
function(mn, un, fudge=100){
   a <- db1(mn,un)
   a[[1]]/(a[[1]] + a[[2]] + fudge)
}
