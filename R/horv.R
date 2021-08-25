# functions for Horvath epigenetic clock
# initial version: 20 Oct 2015 Leo
# Updated version: 18 Jul 2016 TGS 
# Updated once more June 2018 TGS 
# Updated again 12 August 2021 TGS

trafo <- function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }

anti.trafo <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

.split_intercept_from_coeff <- function(x){
  intercept <- min(0, x['(Intercept)'], na.rm = TRUE) 
  interceptless_coeff <- x[names(x) %in% '(Intercept)']
  return(list(intercept=intercept, coeffs=interceptless_coeff))
}

.handle_missing <- function(cpgs, coef_list){
  missing <- names(coef_list$coeffs) %in% names(na.omit(cpgs))
  coef_list$coeffs <- coef_list$coeffs[missing]
  return(coef_list)
}

.calculate_age <- function(x, coeff){
  coef2 <- .split_intercept_from_coeff(coeff)
  coef3 <- .handle_missing(cpgs=x, coef_list=coef2)
  data <- x[names(coef3$coeffs)]
  the_sum <- data %*% coef3$coeffs + coef3$intercept  #data %*% coef2 + coeff[1]
  return(the_sum)
}

.compute_ages <- function(betas, coeff){
  # Calculate on a per-sample basis
  ages <- as.matrix(apply(betas, 2, FUN=.calculate_age, coeff=coeff))
  return(ages)
}
"Horvath"   "Hannum"    "Lin"       "SkinBlood" "PhenoAge" 
agep <- function(betas, coeff = NULL, method = c('horvath', 'hannum', 'phenoage', 'skinblood', 'lin', 'all'), ...){
  data(age_coefficients)
  method <- match.arg(method)
  if(!is.null(coeff)){
    # If coeffs are provided just calculate ages according to provided coefficients
    # Tool is smart enough to find the intercept else if missing will set to 0
    ages <- .compute_ages(betas=betas, coeff=coeff)
  } else {
    ages <- switch(method,
      'horvath' = {
         pre <- .compute_ages(betas=betas, coeff=age_coefficients[['Horvath']])
        # Horvath needs this fancy step
         anti.trafo(pre, adult.age=20)
      },
      'hannum' = {
        .compute_ages(betas=betas, coeff=age_coefficients[['Hannum']])
      },
      'lin' = {
        .compute_ages(betas=betas, coeff=age_coefficients[['Lin']])
      },
      'skinblood' = {
        .compute_ages(betas=betas, coeff=age_coefficients[['SkinBlood']])
      },
      'phenoage' = {
        .compute_ages(betas=betas, coeff=age_coefficients[['PhenoAge']])
      },
      'all' = {
        clocks = c('horvath' = 'Horvath', 'hannum' = 'Hannum', 'phenoage' = 'PhenoAge', 'skinblood' = 'SkinBlood', 'lin' = 'Lin') # Add as many as cases
        do.call('cbind', 
          lapply(clocks, function(x, betas){
            agep(betas = betas, coeff = NULL, method = x)
          }, betas = betas)
        ) # Returns a DF nsample rows, n predictions columns. 
      }
    )
  }
  return(ages)
}

