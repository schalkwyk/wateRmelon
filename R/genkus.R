genkus <-
function( betas, g=getsnp(rownames(betas)) ){ # gets a beta matrix
         apply(betas[g,],1, genkme)
}
