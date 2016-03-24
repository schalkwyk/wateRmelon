\name{outlyx}
\alias{outlyx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify Outliers within Methylumi and Minfi packaged objects
}
\description{
Seeks to identify outliers based on multiple (currently 2) outlier detection methods for (unnormalized) \code{methylumi} and \code{minfi} packaged objects.
}
\usage{
outlyx(x, y = 1, dist1 = 2, dist2 = 0.15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{MethyLumiSet}, \code{MethylSet}, \code{RGChannelSet} object or matrix containing raw betas.
}
  \item{y}{
Desired principal component for outlier identification - only used if other principal components want to be discriminated, only used for IQR outlier detection.
}
  \item{dist1}{
The number of interquartile ranges outliers are to be identified from designated principle component.
}
  \item{dist2}{
Arbitrary cut-off point for identifying outliers via \code{pcout}
}
}
\details{
}
\value{
Returns a character vector containing the names of outliers.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tyler Gorrie-Stone
}
\note{
May perform poorly on normalized data
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 library(wateRmelon)
 data(melon)
 outliers <- outlyx(melon, dist1=2, dist2=0.15)
 # Subsetting Example (won't work with melon dataset)
 # example <- melon[,!colnames(melon)%in%outliers]
}

\keyword{outlier}

