smokp <- function(betas, method='Sst', sst=NULL){
  method <- match.arg(method, choices=c('Sst', 'AHRR', 'McCartney', 'Maas', 'Sugden',
          'Teschendorff', 'Yu', 'Gao', 'Yang', 'Zhang', 'Wen',
          'Langdon',  'Packyears', 'Cessation', 'All'))
  data("smokp_cpgs")
  if(method == 'SSt'){
    intercept <- lapply(smokp_cpgs[c('Never','Former','Current')], function(x){x['(Intercept)']})
    coefs <- lapply(smokp_cpgs[c('Never','Former','Current')], function(x){nonmiss <- x[names(x) %in% rownames(na.omit(betas))]})
    data <- lapply(coefs, function(x){betas[match(names(x), rownames(betas)), , drop = F]})
  } else if(method %in% c('Packyears','Cessation')){
    intercept <- smokp_cpgs[[method]]['(Intercept)']
    coefs <- smokp_cpgs[[method]][names(smokp_cpgs[[method]]) %in% rownames(na.omit(betas))]
    data <- betas[match(names(coefs), rownames(betas)), , drop = F]
  } else if(method == 'Langdon'){
    ilmnids <- split(smokp_cpgs$Langdon[,c('IlmnID','Effect')],smokp_cpgs$Langdon$Classes)
    ilmnids <- lapply(ilmnids, function(x){setNames(x$Effect, x$IlmnID)})
    coefs <- lapply(ilmnids, function(x){nonmiss <- x[names(x) %in% rownames(na.omit(betas))]})
    data <- lapply(coefs, function(x){betas[match(names(x), rownames(betas)), , drop = F]})
  } else {
    coefs <- smokp_cpgs[[method]][names(smokp_cpgs[[method]]) %in% rownames(na.omit(betas))]
    data <- betas[match(names(coefs), rownames(betas)), , drop = F]
  }
  methods <- c('AHRR' = 'ahrr', 'SSt'='sst','Packyears'='smokp','Cessation'='smokp',
               'McCartney'='methylation_score','Maas'='methylation_score',
               'Teschendorff'='smoking_index','Yu'='smoking_index','Gao'='smoking_index','Yang'='smoking_index',
               'Sugden'='sugden','Langdon'='langdon','Wen'='wen','All'='all')
  method_sub <- methods[method]
  smoking <- switch(method_sub,
                    'ahrr' = {
                      out <- data.frame(data['cg05575921',])
                      colnames(out) <- method
                      out
                    },
                    'sst' = {
                      log_odds <- list(
                        Never = as.vector(coefs[['Never']] %*% data[['Never']] + intercept[['Never']]),
                        Former = as.vector(coefs[['Former']] %*% data[['Former']] + intercept[['Former']]),
                        Current = as.vector(coefs[['Current']] %*% data[['Current']] + intercept[['Current']])
                      )
                      odds <- sapply(log_odds, exp)
                      probs <- apply(odds, 2, function(x){x/rowSums(odds)})
                      out <- data.frame(probs, 'SSt' = as.character(factor(max.col(probs), levels=c(1:3), labels=c('Never','Former','Current'))))
                      rownames(out) <- colnames(data$Current)
                      out
                    },
                    'smokp' = {
                      out <- data.frame(t(coefs %*% data + intercept))
                      colnames(out) <- method
                      out
                    },
                    'methylation_score' = {
                      out <- data.frame(t(coefs %*% data)) #in matrix multiplication columns of the left matrix = number of rows of the right
                      colnames(out) <- method
                      out
                    },
                    'smoking_index' = {
                      dataspl <- list(Current = data[,names(sst[sst %in% 'Current'])],
                                      Never = data[,names(sst[sst %in% 'Never'])])
                      means <- lapply(dataspl, function(x){apply(x, 1, mean)})
                      weights <- ifelse((means$Current-means$Never) < 0, -1, 1)
                      sds <- apply(dataspl$Never, 1, sd)
                      diff <- apply(data, 2, function(x){x - means$Never / sds})
                      out <- data.frame(1/length(weights) * colSums(weights * diff))
                      colnames(out) <- method
                      out
                    },
                    'sugden' = {
                      out <- apply(data*100, 2, function(x){x*coefs})
                      out <- data.frame(colMeans(out))
                      colnames(out) <- method
                      out
                    },
                    'langdon' = {
                      evernev <- as.numeric(coefs$'Ever vs never' %*% data$'Ever vs never')
                      curform <- as.numeric(coefs$'Current vs former' %*% data$'Current vs former')
                      out <- cbind(evernev,curform)
                      colnames(out) <- paste0(method, '_', c('ever_never','current_former'))
                      rownames(out) <- colnames(data$'Ever vs never')
                      out
                    },
                    'wen' = {
                      out <- data.frame(1/(1 + exp(-(10.621-10.005*data['cg05575921',] - 8.770*data['cg01940273',]))))
                      colnames(out) <- method
                      out
                    },
                    'all' = {
                      allmethods <- c('SSt','Packyears','Cessation','McCartney','Maas','Teschendorff','Yu','Gao','Yang','Sugden','Langdon','Wen')
                      out <- lapply(allmethods, function(x, betas, sst){
                        smokp(betas = betas, sst = sst, method = x)
                      }, betas=betas, sst=sst)
                      out <- data.frame(do.call('cbind', out))
                    out
                    })
  return(smoking)
}
