\name{agep}
\alias{agep}
\alias{coef}
\title{
Age Prediction from methylomic expression data
}
\description{
Predict age of samples using Horvath's Coefficients
}
\usage{
agep(betas, coeff=NULL, verbose=FALSE)
}
\arguments{
  \item{betas}{
Matrix of betas 
}
  \item{coeff}{
If NULL, will default to Horvath's coefficients. Otherwise can be supplied with 
vector of own coefficients with an intercept as first element.
}
  \item{verbose}{
If TRUE, agep will additionally print out names of probes that are missing from coefficient list within supplied beta matrix.
}
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
\note{
Is extremely sensitive to NAs and will return NA if any probes (listed in coefficients) have NA values.
}
\examples{
data(melon)
agep(betas(melon))
}
