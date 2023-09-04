\name{agep}
\alias{agep}
\alias{agep,MethylSet-method}
\alias{agep,RGChannelSet-method}
\alias{agep,MethyLumiSet-method}
\alias{coef}
\alias{hannumCoef}
\alias{age_coefficients}
\title{
Age Prediction from methylomic expression data
}
\description{
Predict age of samples using Horvaths Coefficients
}
\usage{
agep(betas, coeff = NULL, method = c('horvath', 'hannum', 'phenoage', 'skinblood', 'lin', 'all'), n_missing = TRUE, missing_probes = FALSE, ...)
}
\arguments{
  \item{betas}{
Matrix of betas or MethyLumiSet or MethylSet object.
}
  \item{coeff}{
If NULL, will default to whatever method is specified in method. If not NULL, the expected input should be a vector of coefficients and intercept
}
  \item{method}{
Currently: "horvath", "hannum", "phenoage", "skinblood", "lin" and "all", if "all" agep will seek to calculate ages using all methods else will use the method specified. Default is "horvath".
}
  \item{n_missing}{Logical, additionally output the number of missing CpGs for each sample using the specified method or coeff list.
  }
  \item{missing_probes}{Logical, additionally output the names of missing CpGs for each sample using the specified method or coeff list.
  }
  \item{...}{To pass to arguments to downstream functions to specify adult.age}
}
%\details{}
\value{
Returns matrix of predicted ages per sample. With additional columns created whether n_missing or missing_probes are specified. If method is "all" then all ages will be provided in the same matrix output
}
\references{
Horvath S: DNA methylation age of human tissues and cell types. Genome Biology 2013, 14:R115
}
\author{
Original Functions: Steve Horvath

wateRmelon Implementation: Tyler Gorrie-Stone, Leo Schalkwyk, Louis El Khoury
}
% \note{}

\examples{
data(melon)
agep(melon,coeff=NULL, method="all", n_missing=FALSE)
agep(melon,coeff=NULL, method="horvath", n_missing=TRUE)
}
