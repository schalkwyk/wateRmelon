\name{agep}
\alias{agep}
\alias{agep,MethylSet-method}
\alias{agep,MethyLumiSet-method}
\alias{coef}
\alias{hannumCoef}
\title{
Age Prediction from methylomic expression data
}
\description{
Predict age of samples using Horvath's Coefficients
}
\usage{
agep(betas, coeff=NULL, method=c('horvath','hannum') ,...)
}
\arguments{
  \item{betas}{
Matrix of betas or MethyLumiSet or MethylSet object.
}
  \item{coeff}{
If NULL, will default to Horvath's coefficients. Otherwise can be supplied with 
vector of own coefficients with an intercept as first element.
}
  \item{method}{
currently either 'horvath' or 'hannum'. If using 'hannum' remember to call data(hannumCoef) before analysis and supply coeff with the hannum coefficients, otherwise will default to horvath.
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
