correctII <-
function (BetaValues, SI, SII) 
{
    M <- Beta2M(BetaValues)
    sigma_u <- SII[1]/SI[1]
    sigma_m <- SII[2]/SI[2]
    M <- sapply(M, function(x) {
        if (is.na(x)) return(NA)  ##LS##
        if (x < 0) 
            return(x/sigma_u)
        else return(x/sigma_m)
    })
    return(M2Beta(M))
}
