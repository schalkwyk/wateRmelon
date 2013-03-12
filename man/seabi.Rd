\name{seabi}
\alias{seabi}
\title{
Calculate a performance metric based on male-female differences for Illumina methylation 450K arrays
}
\description{
Calculates an area under ROC curve - based metric for Illumina 450K data using
a t-test for male-female difference as the predictor for X-chromosome location of probes.  The metric is 1-area so that small values indicate good performance, to match our other, standard error based metrics \code{\link{gcose}} and \code{\link{dmrse}}. Note that this requires both male and female samples of known sex and can be slow to compute due to running a t-test on every probe. 
}
\usage{
seabi(bn, stop = 1, sex, X)
}
\arguments{
  \item{bn}{
a matrix of betas (default method) or an object containing betas i.e. a \code{MethyLumiSet} object (\code{methylumi} package), a \code{MethylSet} or \code{RGChannelSet} object (\code{minfi} package) or a \code{exprmethy450} object (\code{IMA} package).

}
  \item{stop}{
partial area under curve is calculated if stop value <1 is provided
}
  \item{sex}{
a factor giving the sex of each sample (column)
}
  \item{X}{
a logical vector of length equal to the number of probes, true for features mapped to X-chromosome 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{

a value between 0 and 1.  values close to zero indicate high data quality as judged by the ability to discriminate male from female X-chromosome DNA methylation.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
}
\author{
leonard.schalkwyk@kcl.ac.uk
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
   library(methylumi)
   data(melon)
   sex  <- pData(melon)$sex
   X    <- melon@featureData@data$CHR=='X'
   seabi(betas(melon), sex=sex, X=X)

# methylumi method
   seabi(melon, sex=sex, X=X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
