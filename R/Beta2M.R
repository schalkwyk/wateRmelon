#' Internal functions for peak.correction (fuks)
#' 
#' Internal functions for peak.correction
#' 
#' 
#' @aliases Beta2M M2Beta correctI correctII summits
#' @param B %% ~~Describe \code{B} here~~ a vector or matrix of beta values for
#' conversion
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... a vector or matrix of the same shape as the input
#' @author Matthieu Defrance <defrance@@bigre.ulb.ac.be>
#' @references
#' 
#' Dedeurwaerder S, Defrance M, Calonne E, Sotiriou C, Fuks F: Evaluation of
#' the Infinium Methylation 450K technology . Epigenetics 2011, 3(6):771-784.
#' @export Beta2M
Beta2M <-
function (B) 
{
    return(log2(B/(1 - B)))
}
