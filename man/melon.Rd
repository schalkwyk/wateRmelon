\name{melon}
\alias{melon}
\docType{data}
\title{
Small MethyLumi data set for examples and testing
}
\description{
This object was derived using \code{\link[methylumi]{methylumiR}} on an 
edited GenomeStudio file containing a small subset of features.
It works with all of the wateRmelon package beta functions
(see \code{\link{dasen}} and metrics (see \code{\link{genki}}, \code{\link{seabi}}, and \code{\link{dmrse_col}}) except for \code{\link{swan}}.
}
\usage{data(melon)}
\format{
MethyLumiSet with assayData containing 3363 features, 12 samples 
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
}
%\references{
%}
\examples{
library(methylumi)
data(melon)
boxplot(log(methylated(melon)), las=2)
## maybe str(melon) ; plot(melon) ...
}
\keyword{datasets}
\keyword{MethyLumiSet}
