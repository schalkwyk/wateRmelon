dfs2 <-
function(x, onetwo){  # x is a matrix of intensities
                              # onetwo is a character vector 
                              # of same order and length 
                              # indicating assay I or II 
   one <- density(x[onetwo=='I'], na.rm=T, n = 2^15, from = 0, to = 5000)
   two <- density(x[onetwo=='II'],na.rm=T, n = 2^15, from = 0, to = 5000)
   one$x[which.max(one$y)] - two$x[which.max(two$y)]
}
