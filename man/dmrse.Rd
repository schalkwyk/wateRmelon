\name{dmrse}
\alias{dmrse}
\alias{dmrse_col}
\alias{dmrse_row}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard error of iDMR 450k array DNA methylation features
}
\description{Imprinting differentially methylated regions (iDMRs) are expected to be approximately half methylated, as is observed at the 227 probes in known iDMRs. These functions calculate measures of dispersion for the beta values at these CpG sites, of which the most useful is dmrse_row, which is the between-sample 
standard error.
}
\usage{
dmrse(betas, idmr = iDMR())
dmrse_col(betas, idmr = iDMR())
dmrse_row(betas, idmr = iDMR())
}
\arguments{
  \item{betas}{
a matrix of betas (default method), a \code{MethyLumiSet} object (\code{methylumi} package), a \code{MethylSet} or \code{RGChannelSet} object (
\code{minfi} package) or a \code{exprmethy450} object (\code{IMA} package).
}
  \item{idmr}{
%%     ~~Describe \code{idmr} here~~
a character vector of iDMR probe names such as returned by iDMR() 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
return a standard error of the mean of betas for all samples and 
iDMR probes (dmrse) or the standard error of the mean for just the between sample component(dmrse_row) or between probe(dmrse_col) component.
}
\references{
Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
}
\author{
%%  ~~who you are~~
Leonard.Schalkwyk@kcl.ac.uk
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{seabi}}, a sex-difference metric, and 
\code{\link{genki}}, based on SNPs.

}
\examples{

  #MethyLumiSet method
     data(melon)
     dmrse(melon)

  #MethyLumiSet method after normalization
     melon.dasen <- dasen(melon)
     dmrse(melon.dasen)
}

