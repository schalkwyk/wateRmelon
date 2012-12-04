\name{metrics}
\alias{metrics}
\title{
Calculate a full set of 450K normalization/performance metrics
}
\description{
Calculate X-chromosome, SNP and imprinting DMR metrics for a matrix of betas
from an Illumina 450K Human DNA methylation array.  Requires precalculated
t-test p-values for sex differences, a list of X-chromosome features and 
of imprinting DMR features.}
\usage{
metrics(betas, pv, X, idmr = iDMR, subset = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{
a matrix of betas, each row representing a probe, each column a sample
}
  \item{pv}{
a vector of p-values such as produced by \code{sextest}, one per row of betas
}
  \item{X}{
a logical vector of the same length as \code{pv}, indicating whether each probe
is mapped to the X-chromosome
}
  \item{idmr}{
a character vector of probe names known to be in imprinting DMRs.  Can be obtained with \code{iDMR()} or \code{data(iDMR)}
}
  \item{subset}{
index or character vector giving a subset of betas to be tested
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{ dmrse_row}{see \code{dmrse_row}}
\item{ dmrse_col}{see \code{dmrse_col}}
\item{ dmrse    }{see \code{dmrse    }}
\item{ gcoms_a  }{see \code{genki    }}
\item{ gcose_a  }{see \code{genki    }}
\item{ gcoms_b  }{see \code{genki    }}
\item{ gcose_b  }{see \code{genki    }}
\item{ gcoms_c  }{see \code{genki    }}
\item{ gcose_c  }{see \code{genki    }}
\item{ seabird  }{see \code{seabi    }}
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
%%  ~~further notes~seabi~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(melon)
melon.dasen <- dasen(melon)
bn <-betas(melon.dasen)
X <- melon.dasen@featureData@data$CHR=='X'
data(iDMR)
sex <- pData(melon.dasen)$sex
pv <- sextest(bn,sex)
melon.metrics <- metrics(bn, pv, X, idmr = iDMR, subset = NULL) 
}
