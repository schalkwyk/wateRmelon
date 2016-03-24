
# functions for Horvath epigenetic clock
# initial version 20 Oct 2015 Leo
# Made clone version: 16 Mar 2016 TGS 'agep2'
#   re: Handling Missing Data + allow for custom coefficients. 

trafo <- function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

agep <- function(betas){
 data(coef)
 data <- betas[names(coef)[-1],]
# some kind of error handling re missing data should go here
# horvath imputes using his 'gold standard'
   pre <- t(data) %*% coef[-1] + coef[1]
   anti.trafo(pre)
}

# Use if one wishes to use own coefficients.
agep2 <- function(betas, y = NULL){
 if(is.null(y)){data(coef)
                y <- coef}
 data <- betas[names(y)[-1],]
 pre <- t(data) %*% y[-1] + y[1]
 anti.trafo(pre)
}

# Attempt to account for missing values - currently not working!!
agep3 <- function(betas, y = NULL, verbose=F){
# If one wishes to use their own set of coefficients.
# Only way I could figure out how to call data(coef) without
# throwing up an error
 if(is.null(y)){data(coef)
                y <- coef}
# This is one way to avoid NAs in the output albeit 'slow'.
# Determining missing probes
# betas <- na.omit(betas)
 miss <- rownames(betas)[rownames(betas)%in%names(y)]
 coef2 <- y[names(y)%in%miss]
 data <- betas[miss,]
# Per louis' request:
 if(verbose){cat(((length(y)-1)-length(miss)), 'probes missing:','\n')
             if(length(miss)>0){ print(names(y)[-1][!names(y)[-1]%in%miss])}
             cat('\n')}
 pre <- t(data) %*% coef2 + y[1]
 anti.trafo(pre)
}


