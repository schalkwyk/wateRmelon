sextest <-
function(betas, sex, ...) {# ... might be useful for subsets
   mod2 <- function(x){ 
      r <- try(t.test(x ~ sex, ... ), silent=TRUE)
      if(inherits(r,'try-error')) return(NA)
      r$p.
   }
   apply(betas,1,mod2, ...)
}
