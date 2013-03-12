seabird <-
function(pr, stop=1, X){ # sexdiff pvalue ROC AUC
                                    # pr : pvals from sextest
                                    # stop: fraction for partial AUC
                                    # X: logical vector-probe on X?
# ROCR version
#    pr <- prediction(1 - pr, X)
#    unlist(performance(pr, "auc", fpr.stop = stop)@y.values)
# pROC version
#    as.numeric(auc(formula=X~pr, data=NULL)) 
# ROC version (bioconductor)

   pAUC(rocdemo.sca(truth=X, data=1-pr), stop) 

}
   



seabi <- 
function (bn, stop=1, sex, X){
   1 - seabird( sextest(bn, sex), stop, X )
}
