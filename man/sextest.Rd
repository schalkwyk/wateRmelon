\name{sextest}
\alias{sextest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test Illumina methylation 450K array probes for sex difference (internal function for calculating
\code{seabi} performance metric)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a wrapper for \code{lm} which does the equivalent of a Student t-test for difference in betas between males and females for each row of a matrix of betas.
}
\usage{
sextest(betas, sex, ...)
}
\arguments{
  \item{betas}{
a matrix of betas, each row is a probe, each column a sample
}
  \item{sex}{
a factor with 2 levels for male and female
}
  \item{\dots}{
additional arguments to be passed to \code{lm}  
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a vector of p-values of length equal to the number of rows of betas
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
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{seabi}}
\code{\link{seabird}}
}
\examples{

   #MethyLumiSet method
     data(melon)
     sex  <- pData(melon)$sex
     melon.sextest<-sextest(betas(melon),sex)


   #MethyLumiSet method with quality control step
     data(melon)
     melon.dasen <- dasen(melon)
     sex  <- pData(melon.dasen)$sex
     melon.sextest<-sextest(betas(melon.dasen),sex)
}
