\name{pwod}
\alias{pwod}
\title{
Probe-Wise Outlier Detection
}
\description{
'P'robe-'W'ise 'O'utlier 'D'etection via interquartile ranges.
 }
\usage{
pwod(object, mul=4)
}
\arguments{
  \item{object}{MethyLumiSet, RGChannelSet, MethylSet object or matrix containing betas.
}
  \item{mul}{
Number of interquartile ranges used to determine outlying probes. Default is 4 to ensure only very
obvious outliers are removed.
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
new.betas <- pwod(cattle, mul=4)
}




