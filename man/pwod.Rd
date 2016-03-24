\name{pwod}
\alias{pwod}
\title{
Probe-Wise Outlier Detection
}
\description{
'P'robe-'W'ise 'O'utlier 'D'etection via interquartile ranges.
 }
\usage{
library(wateRmelon)
data(melon)
cattle <- betas(melon)
pwod(cattle) -> filter 
}
\arguments{
  \item{betas}{
Matrix of betas 
}
  \item{mul}{
Number of interquatile ranges used to determine outlying probes.
}
}
\details{
Detects outlying probes across arrays in methylumi and minfi objects. Outliers are probable low MAF SNP heterozygotes(?)
}

\value{
Returns inputted beta matrix with outlying probes coerced to NA
}
\references{

}

\author{
Tyler Gorrie-Stone - tgorri@essex.ac.uk
}

\note{

}

