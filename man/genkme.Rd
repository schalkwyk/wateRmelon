\name{genkme}
\alias{genkme}
\alias{genkus}
\alias{getsnp}
\alias{gcose}
\alias{gcoms}
\title{
Internal functions for genotype-based normalization metrics
}
\description{
genkme - genotype calling with 1d k-means

genkus - apply genkme to available SNPs

getsnp - grep the rs-numbered probes

gcose - calculate between-sample SNP standard error

gcoms - calculate between-sample SNP mean-squared deviation
}
\usage{
genkme(y, peaks = c(0.2, 0.5, 0.8))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a vector or matrix of numeric values (betas, between 0 and 1)
}
  \item{peaks}{
initial values for cluster positions
}
}
\details{
see \code{\link{genki}}
}
\value{

see \code{\link{genki}}
}
\references{
Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
}
\author{
Leonard.Schalkwyk@kcl.ac.uk
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

