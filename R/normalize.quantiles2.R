normalize.quantiles2 <-
function(X, Reference.Quantiles){

    apply(X, 2, function(x, y) y[rank(x)], Reference.Quantiles)
}
