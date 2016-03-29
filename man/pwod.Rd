\name{pwod}
\alias{pwod}
\title{
Probe-Wise Outlier Detection
}
\description{
'P'robe-'W'ise 'O'utlier 'D'etection via interquartile ranges.
 }
\usage{
pwod(object, mul)
}
\arguments{
  \item{object}{MethyLumiSet, RGChannelSet, MethylSet object or matrix containing betas.
}
  \item{mul}{
Number of interquatile ranges used to determine outlying probes.
}

}
\details{
Detects outlying probes across arrays in methylumi and minfi objects. Outliers are probable low MAF/SNP heterozygotes.
}
\value{
Returns supplied beta matrix with outlying probes coerced to NA
}
%\references{}
%\note{}
\author{
Tyler Gorrie-Stone - tgorri@essex.ac.uk
}
\examples{
library(wateRmelon)
data(melon)
cattle <- betas(melon)
newbetas <- pwod(cattle, mul=4)
}




