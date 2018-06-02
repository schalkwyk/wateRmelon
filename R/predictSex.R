predictSex <- function(x, x.probes, pc = 2, plot=FALSE, irlba=TRUE, center=FALSE, scale. = FALSE){
  x <- x[x.probes,]
  cm <- colMedians(x) # colMeans?
  message('Performing PCA')
  if(irlba & require('irlba')) pca <- prcomp_irlba(x, n = pc, center = center, retx = F, scale. = scale.)$rot 
  else pca <- prcomp(x, center = center, retx = F, scale. = scale.)$rot[,seq_len(pc)]
  if(any(table(cm >= .48) >= (ncol(x)*.95))){
    message('Data likely consists of single gender')
  }   
  message('Generating Clusters')
  c1 <- kmeans(cm, centers=fivenum(cm)[c(2,4)])
  c2 <- kmeans(pca[,pc], centers=2)
  n <- 0
  while(sum((s <- c1$cluster+c2$cluster)==3) > (length(c2$cluster)/2)){
    c2 <- kmeans(pca[,pc], centers=2)
    n <- n + 1
    if(n > 10) stop('Cannot find meaningful clusters for data')
  }
  is.na(s[s==3]) <- TRUE  
  o <- order(t <- tapply(cm, s, mean))
  s[s==names(t)[o[1]]] <- 'Male'
  s[s==names(t)[o[2]]] <- 'Female'
  s[is.na(s)] <- 'Undefined'
  if(plot){
    par(mar=c(4.5,4.5,1,8))
    plot(x=pca[,1], y=pca[,pc], col=factor(s), xlab="PC 1", ylab=paste('PC', pc))
    legend('right', legend=levels(factor(s)), col=1:3, pch=1,inset=c(-0.3,0),xpd=T) 
  }
  return(s)
}
