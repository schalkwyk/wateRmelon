smokp <- function(betas, method, sst=NULL){
  load('smokp_cpgs.Rda')
  if(method == 'SSt'){
    intercept <- lapply(smokp_cpgs[c('Never','Former','Current')], function(x){x['(Intercept)']})
    coefs <- lapply(smokp_cpgs[c('Never','Former','Current')], function(x){
      miss <- names(x) %in% rownames(na.omit(betas))
      t(as.matrix(x[miss]))})
    data <- lapply(coefs, function(x){betas[match(colnames(x), rownames(betas)),]})
  } else {
    intercept <- smokp_cpgs[[method]]['(Intercept)']
    miss <- names(smokp_cpgs[[method]]) %in% rownames(na.omit(betas))
    coefs <- t(as.matrix(smokp_cpgs[[method]][miss]))
    data <- betas[match(colnames(coefs), rownames(betas)),]
  }
  methods <- data.frame(method_sub=c('SSt',rep('smokp',2), rep('Methylation_Score',3),rep('Smoking_Index',4),'Zhang','Sugden'),
                        row.names=c('SSt','Packyears','Cessation','McCartney','Christiansen','Odintsova',
                                    'Teschendorff','Yu','Gao','Yang','Zhang','Sugden')
                        , stringsAsFactors = F)
  method_sub <- methods[method,'method_sub']
  out <- switch(method_sub,
                'SSt' = {
                  log_odds <- list()
                  for(i in c('Never','Former','Current')){
                    log_odds[[i]] <- as.vector(coefs[[i]] %*% data[[i]] + intercept[[i]])
                  } 
                  odds <- sapply(log_odds, exp)
                  log_odds <- data.frame(do.call(cbind,log_odds))
                  sum_odds <- rowSums(odds)
                  probs <- data.frame(apply(odds, 2, function(x){x/sum_odds}))
                  out <- data.frame(as.character(factor(max.col(probs), levels=c(1:3), labels=c('Never','Former','Current'))))
                },
                'smokp' = {
                  out <- data.frame(t(coefs %*% data + intercept))
                },
                'Methylation_Score' = {
                  out <- data.frame(t(coefs %*% data))
                },
                'Smoking_Index' = {
                  nevdat <- data[,match(names(sst[sst %in% 'Never']), colnames(data))]
                  curdat <- data[,match(names(sst[sst %in% 'Current']), colnames(data))]
                  nev_means <- apply(nevdat, 1, mean)
                  cur_means <- apply(curdat, 1, mean)
                  weights <- cur_means-nev_means
                  weights <- ifelse(weights < 0, -1, 1)
                  weights <- weights[colnames(coefs)]
                  sds <- apply(nevdat, 1, sd)
                  out <- apply(data, 2, function(x){x - nev_means / sds})
                  out <- colSums(weights * out)
                  out <- data.frame(1/length(weights) * out)
                },
                'Zhang' = {
                  data <- t(as.matrix(data)) #only 1 Zhang CpG in beta matrix
                  rownames(data) <- colnames(coefs)
                  lowest_quartile <- setNames(apply(data,1,quantile)[2,],colnames(coefs))
                  out <- list()
                  for(i in colnames(coefs)){
                    out[[i]] <- ifelse(data[i,] < lowest_quartile[i], 1, 0)
                  }
                  out <- data.frame(Reduce(`+`, out))
                },
                'Sugden' = {
                  data <- data*100
                  out <- apply(data, 2, function(x){x*coefs})
                  out <- data.frame(colMeans(out))
                })
  colnames(out) <- method
  return(out)
}
