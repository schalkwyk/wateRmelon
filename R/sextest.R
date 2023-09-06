#' Test Illumina methylation 450K array probes for sex difference (internal
#' function for calculating \code{seabi} performance metric)
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ This
#' is a wrapper for \code{lm} which does the equivalent of a Student t-test for
#' difference in betas between males and females for each row of a matrix of
#' betas.
#' 
#' 
#' @param betas a matrix of betas, each row is a probe, each column a sample
#' @param sex a factor with 2 levels for male and female
#' @param \dots additional arguments to be passed to \code{lm}
#' @return Returns a vector of p-values of length equal to the number of rows
#' of betas
#' @author lschal@@essex.ac.uk
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link{seabi}} \code{\link{seabird}}
#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted)
#' @examples
#' 
#' 
#'    #MethyLumiSet method
#'      data(melon)
#'      sex  <- pData(melon)$sex
#'      melon.sextest<-sextest(betas(melon),sex)
#' 
#' 
#'    #MethyLumiSet method with quality control step
#'      data(melon)
#'      melon.dasen <- dasen(melon)
#'      sex  <- pData(melon.dasen)$sex
#'      melon.sextest<-sextest(betas(melon.dasen),sex)
#' 
#' @export sextest
sextest <-
function(betas, sex, ...) {# ... might be useful for subsets
   mod2 <- function(x){ 
      r <- try(t.test(x ~ sex, ... ), silent=TRUE)
      if(inherits(r,'try-error')) return(NA)
      r$p.
   }
   apply(betas,1,mod2, ...)
}
