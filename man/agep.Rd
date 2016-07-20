\name{agep}
\alias{agep}
\alias{agep,MethylSet-method}
\alias{agep,MethyLumiSet-method}
\alias{coef}
\title{
Age Prediction from methylomic expression data
}
\description{
Predict age of samples using Horvath's Coefficients
}
\usage{
agep(betas, coeff=NULL, verbose=FALSE,...)
}
\arguments{
  \item{betas}{
Matrix of betas or MethyLumiSet or MethylSet object.
}
  \item{coeff}{
If NULL, will default to Horvath's coefficients. Otherwise can be supplied with 
vector of own coefficients with an intercept as first element.
}
  \item{verbose}{
If TRUE, agep will additionally print out names of probes per sample that are missing from coefficient list within supplied beta matrix.
}
  \item{...}{To pass to arguments to downstream functions to specify adult.age}
}
%\details{}
\value{
Returns matrix of predicted ages per sample.
}
\references{
Horvath S: DNA methylation age of human tissues and cell types. Genome Biology 2013, 14:R115
}
\author{
Original Functions: Steve Horvath

wateRmelon Implementation: Leo Schalkwyk, Louis El Khoury
}
% \note{}

\examples{
data(melon)
agep(melon,coeff=NULL)
# or
agep(betas(melon),coeff=NULL)
}
