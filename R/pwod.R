 # Author: Tyler Gorrie-Stone, tgorri@essex.ac.uk
 # Revision Date: 07-01-2016
pwod <- function(object, mul=4){
 # 'P'robe-'W'ise 'O'utlier 'D'etection via interquartile ranges.
 # -- probable low MAF SNP heterozygotes --
 # Arguments:
 #  object   : methylumi object or minfi object or Matrix of Betas
 #  mul     : Number of interquartile ranges to determine outlying probes.
 #            Default = 4. To exclude the very obvious.
 #
 # Returns  : Matrix of Betas with outlying probes coerced to NA.
  filter <- t(apply(object, 1, function(x){
    quan <- fivenum(x)
    iqr <- quan[4] - quan[2]
    bounds <- c(quan[4] + (iqr * mul),
                quan[2] - (iqr * mul)) 
  # Upper Bound is [1], Lower Bound is [2]
    d <- x > bounds[1] | x < bounds[2]
    x[d] <- NA 
    return(x)
    }))
 # Calculating total number of outlying probes.
 # To give an idea of what has changed (if anything).
  tot <- sum(is.na(filter)) - sum(is.na(object))
  cat(tot,"probes detected.", "\n") 
 # Output
  filter
}


