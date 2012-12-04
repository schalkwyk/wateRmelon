daten2 <-
function(mn, un, onetwo, fudge=100, ...){

  mnc <- dfsfit(mn, onetwo, ...)
  unc <- dfsfit(un, onetwo, ...)
  normalizeQuantiles(mnc) -> mncn
  normalizeQuantiles(unc) -> uncn
  mncn/(mncn + uncn + fudge)

}
