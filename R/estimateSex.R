#' Predict sex by using robust sex-related CpG sites on ChrX and ChrY
#'
#' @param betas raw beta values, ideally: beta = M / (M + U + 100).
#' @param do_plot logical. Should plot the predicted results? Default: FALSE
#'
#' @return dataframe contains predicted sex information.
#' @export
#' @author Wang, Yucheng
#'
#' @examples
#' pred_XY <- estimateSex(betas, do_plot=TRUE)
estimateSex <- function(betas, do_plot=FALSE){
  # predict sex by two PCAs on X and Y chromosomes
  data("sexCoef")
  # Z score normalization
  betas <- betas[rownames(betas) %in% sex_coef$IlmnID, ]
  message('Normalize beta values by Z score...')
  autosomes <- sex_coef$IlmnID[!(sex_coef$CHR %in% c('X', 'Y'))]
  auto_betas <- betas[rownames(betas) %in% autosomes, ]
  d_mean <- colMeans(auto_betas, na.rm=TRUE)
  d_sd <- colSds(auto_betas, na.rm=TRUE)
  z_beta <- (t(betas) - d_mean) / d_sd
  message('Fishished Zscore normalization.')

  # Sex prediction
  pred_XY <- list()
  for(chr in c('X', 'Y')){
    coefs <- sex_coef[sex_coef$pca == chr,]
    chr_beta <- z_beta[, coefs$IlmnID]
    chr_beta[is.na(chr_beta)] <- 0
    pred_chr <- t(t(chr_beta) - coefs$mean) %*% coefs$coeff
    pred_XY[[chr]] <- pred_chr
  }
  pred_XY <- data.frame(pred_XY)

  pred_XY$'predicted_sex' <- 'Female'
  pred_XY$'predicted_sex'[(pred_XY$X < 0) & (pred_XY$Y > 0)] <- 'Male'
  pred_XY$'predicted_sex'[(pred_XY$X > 0) & (pred_XY$Y > 0)] <- '47,XXY'
  pred_XY$'predicted_sex'[(pred_XY$X < 0) & (pred_XY$Y < 0)] <- '45,XO'
  if(do_plot){
    plot_predicted_sex(pred_XY)
  }else{
    message('You can visualize the predicted results by set "do_plot=TRUE".\n')
  }
  return(pred_XY)
}


plot_predicted_sex <- function(pred_XY){
  # visualization of predicted sex
  plot(Y~X, data=pred_XY, pch=1, xlab='ChrX-PC1', ylab='ChrY-PC1')
  abline(v=0, lty='dashed')
  abline(h=0, lty='dashed')
  abnormls <- pred_XY[!(pred_XY$'predicted_sex' %in% c('Male', 'Female')),]
  if(nrow(abnormls) > 0){
    points(Y~X, data=abnormls, pch=2, col='red')
    for(i in 1:nrow(abnormls)){
      text(abnormls$X[i], abnormls$Y[i], rownames(abnormls)[i], pos=3, col='red', cex=0.5)
    }
  }
  text(-10, 2, '46,XY', cex=1.2, col='blue')
  text(-10, -2, '45,XO', cex=1.2, col='blue')
  text(10, -2, '46,XX', cex=1.2, col='blue')
  text(10, 2, '47,XXY', cex=1.2, col='blue')

}
