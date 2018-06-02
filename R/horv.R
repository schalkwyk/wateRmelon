# functions for Horvath epigenetic clock
# initial version: 20 Oct 2015 Leo
# Updated version: 18 Jul 2016 TGS 
#    Now Handles missing probes based on single columns.

trafo <- function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

agep <- function(betas, coeff = NULL, method = c('horvath', 'hannum')){
  method <- match.arg(method)
  if(method == 'hannum' & is.null(coeff)) stop('Please supply coeffs for hannum\'s clock.')
  if(is.null(coeff)){
    message('No coefficients detected, using Horvaths clock')
    data(coef)
    coeff <- coef
  }
  if(method == 'horvath'){
    ages <- as.matrix(apply(betas,2,function(x){
      miss <- names(coeff)[-1]%in%names(na.omit(x))
      coef2 <- coeff[-1][miss]
      data <- x[names(coef2)]
      pre <- data %*% coef2 + coeff[1]
      anti.trafo(pre, adult.age=20) 
    }))
  } else if(method == 'hannum'){
    ages <- as.matrix(apply(betas,2,function(x){
      miss <- names(coeff)%in%names(na.omit(x))
      coef2 <- coeff[miss]
      data <- x[names(coef2)]
      data %*% coef2 + 0 #?  
    }))
  }
  return(ages)
}

