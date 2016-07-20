# functions for Horvath epigenetic clock
# initial version: 20 Oct 2015 Leo
# Updated version: 18 Jul 2016 TGS 
#    Now Handles missing probes based on single columns.

trafo <- function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

agep <- function(betas, coeff = NULL, verbose = FALSE,...){
  # If no coeff is supplied:
  if(is.null(coeff)){
    data(coef)
    coeff <- coef
  }

  as.matrix(apply(betas, 2, function(x){
    miss <- names(coeff)[-1]%in%names(na.omit(x)) # Identify missing probes in sample
    coef2 <- coeff[-1][miss] # Use only probes that are present
    # Subset Betas by Probes
    data <- x[names(coef2)]
    # Per louis request: Prints out missing probes...
      if(verbose){
        if(sum(!miss)>0){
          cat(sum(!miss), 'probes missing:','\n')
          print(names(coeff)[-1][!miss]) # Ugly output - consider changing. 
          cat('\n')
        }
      }
    # Age prediction
    pre <- data %*% coef2 + coeff[1]
    anti.trafo(pre, adult.age=20) 
    }
    )
  ) 
} # OK?
