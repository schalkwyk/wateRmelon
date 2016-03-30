\name{outlyx}
\alias{outlyx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify Outliers within Methylumi and Minfi packaged objects
}
\description{
Seeks to identify outliers based on multiple (currently 2) outlier detection methods for \code{methylumi} and \code{minfi} packaged objects.
}
\usage{
outlyx(x, iqr=TRUE, iqrP=2, pc=1,
       mv=TRUE, mvP=0.15, plot=FALSE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{MethyLumiSet}, \code{MethylSet}, \code{RGChannelset} object or matrix containing raw betas.
}
  \item{iqr}{
If TRUE, the outliers based on interquartile ranges will be determined
}
  \item{iqrP}{
The number of interquartile ranges outliers are to be identified from designated principle component.
}
  \item{pc}{
Desired principal component for outlier identification - only used if other principal components want to be discriminated, only used for IQR outlier detection.
}
  \item{mv}{
If TRUE, the outliers will detected using \code{pcout}
}
  \item{mvP}{
Arbitrary cut-off point for identifying outliers via \code{pcout}
}
  \item{plot}{
If TRUE, alongside regular output, a plot will be constructed displaying relative 'location'
of each sample. Outliers are those that fall within the highlighted regions.
}
  \item{...}{
Additional arguments passed to \code{pcout}
}
}
%\details{}
\value{
Returns a dataframe of TRUE/FALSE per sample. Where TRUE is outlying.
}
%\references{}
\author{
Tyler Gorrie-Stone - tgorri@essex.ac.uk
}
\note{
May perform poorly on normalized data
}
\examples{
 library(wateRmelon)
 data(melon)
 outliers <- outlyx(melon,iqr=TRUE, iqrP=2, pc=1,
                    mv=TRUE, mvP=0.15, plot=TRUE)
}

\keyword{outlier}
