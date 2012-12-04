naten <-
function (mn, un, fudge=100 ) {

   mnn <- normalizeQuantiles(mn)
   unn <- normalizeQuantiles(un)
   mnn/(mnn + unn + fudge)

}
