seabird <-
function(pr, stop=1, X){ # sexdiff pvalue ROC AUC
                                    # pr : pvals from sextest
                                    # stop: fraction for partial AUC
                                    # X: logical vector-probe on X?
   if(!library(ROCR, logical.return=TRUE, quietly=TRUE)){
      stop('can\'t load ROCR package')
   }
   pr   <- prediction(1-pr , X)
   unlist(performance(pr, "auc", fpr.stop=stop)@y.values)
}


seabi <- 
function (bn, stop=1, sex, X){
   1 - seabird( sextest(bn, sex), stop, X )
}
