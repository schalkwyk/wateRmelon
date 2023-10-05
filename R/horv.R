# functions for Horvath epigenetic clock
# initial version: 20 Oct 2015 Leo
# Updated version: 18 Jul 2016 TGS 
# Updated once more June 2018 TGS 
# Updated again 12 August 2021 TGS

trafo <- function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }

anti.trafo <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

.split_intercept_from_coeff <- function(x){
  intercept <- x['(Intercept)']
  if(is.na(intercept)) intercept <- 0
  interceptless_coeff <- x[!names(x) %in% '(Intercept)']
  return(list(intercept=intercept, coeffs=interceptless_coeff))
}

.handle_missing <- function(cpgs, coef_list){
  not_missing <- names(coef_list$coeffs) %in% names(na.omit(cpgs))
  coef_list$missing_probes <- paste0(na.omit(names(coef_list$coeffs)[!not_missing]), collapse=';')
  coef_list$n_missing <- length(na.omit(names(coef_list$coeffs)[!not_missing]))
  coef_list$coeffs <- coef_list$coeffs[not_missing]
  return(coef_list)
}

.calculate_age <- function(x, coeff){
  coef2 <- .split_intercept_from_coeff(coeff)
  coef3 <- .handle_missing(cpgs=x, coef_list=coef2)
  data <- x[names(coef3$coeffs)]
  the_sum <- data %*% coef3$coeffs + coef3$intercept  #data %*% coef2 + coeff[1]
  return( data.frame(the_sum, coef3$n_missing, coef3$missing_probes, stringsAsFactors = FALSE) )
}

.compute_ages <- function(betas, coeff){
  # Calculate on a per-sample basis
  ages <- apply(betas, 2, FUN=.calculate_age, coeff=coeff)
  ages <- do.call('rbind', ages)
  return(ages)
}

agep <- function(betas, coeff = NULL, method = c('horvath', 'hannum', 'phenoage', 'skinblood', 'lin', 'all'), n_missing = TRUE, missing_probes = FALSE, ...){
  data("age_coefficients")
  method <- match.arg(method)
  if(!is.null(coeff)){
    # If coeffs are provided just calculate ages according to provided coefficients
    # Tool is smart enough to find the intercept else if missing will set to 0
    ages <- .compute_ages(betas=betas, coeff=coeff)
    colnames(ages) <- c('custom_age', 'n_missing', 'missing_probes')
  } else {
    ages <- switch(method,
      'horvath' = {
         pre <- .compute_ages(betas=betas, coeff=ageCoefs[['Horvath']])
         # Horvath needs this fancy step
         pre[,1] <- anti.trafo(pre[,1], adult.age=20)
         colnames(pre) <- c('horvath.age', 'horvath.n_missing', 'horvath.missing_probes')
         pre
      },
      # Prime the rest in switches incase we need to do more things to each individually...
      'hannum' = {
        pre <- .compute_ages(betas=betas, coeff=ageCoefs[['Hannum']])
        colnames(pre) <- c('hannum.age', 'hannum.n_missing', 'hannum.missing_probes')
        pre
      },
      'lin' = {
        pre <-.compute_ages(betas=betas, coeff=ageCoefs[['Lin']])
        colnames(pre) <- c('lin.age', 'lin.n_missing', 'lin.missing_probes')
        pre
      },
      'skinblood' = {
        pre <- .compute_ages(betas=betas, coeff=ageCoefs[['SkinBlood']])
        pre[,1] <- anti.trafo(pre[,1], adult.age=20) # I think...
        colnames(pre) <- c('skinblood.age', 'skinblood.n_missing', 'skinblood.missing_probes')
        pre
      },
      'phenoage' = {
        pre <- .compute_ages(betas=betas, coeff=ageCoefs[['PhenoAge']])
        colnames(pre) <- c('phenoage.age', 'phenoage.n_missing', 'phenoage.missing_probes')
        pre
      },
      'all' = {
        clocks = c('horvath' = 'horvath', 'hannum' = 'hannum', 'phenoage' = 'phenoage', 'skinblood' = 'skinblood', 'lin' = 'lin') # Add as many as cases
        out <-
          lapply(clocks, function(x, betas){
            agep(betas = betas, coeff = NULL, method = x, n_missing=n_missing, missing_probes=missing_probes )
          }, betas = betas)
        out <- do.call('cbind', out)
        return(out)
      }
    )
  }
  return(ages[,c(TRUE, n_missing, missing_probes)])
}


