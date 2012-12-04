summits <-
function (BetaValues)
{
    d <- density(BetaValues, na.rm=T)  # LS na.rm
    yneg <- d$y[1:which(d$x > M2Beta(0))[1]]
    ypos <- d$y[which(d$x > M2Beta(0))[1]:length(d$y)]
    sa <- d$x[which(d$y == max(yneg))]
    sb <- d$x[which(d$y == max(ypos))]
    return(c(Beta2M(sa), Beta2M(sb)))
}
