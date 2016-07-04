# functions for Horvath epigenetic clock
# initial version: 20 Oct 2015 Leo
# Updated version: 16 Mar 2016 TGS 
#   re: Handling Missing Data + allow for custom coefficients. 

trafo <- function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

agep <- function(betas, coeff = NULL, verbose=FALSE){
# If no coeff is supplied:
 if(is.null(coeff)){
   data(coef)
   coeff <- coef
   }
# Detect Missing Probes from beta matrix
 miss <- names(coeff)[-1]%in%rownames(betas)
 coef2 <- coeff[-1][miss]
# Subset Betas by Probes
 data <- betas[names(coef2),]
# Per louis request: Prints out missing probes
 if(verbose){
   if(sum(!miss)>0){
     cat(sum(!miss), 'probes missing:','\n')
     print(names(coeff)[-1][!miss])
     cat('\n')
     }
   }
# Age prediction
 pre <- t(data) %*% coef2 + coeff[1]
 anti.trafo(pre)
}

