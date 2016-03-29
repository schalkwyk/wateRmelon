\name{qual}
\alias{qual}
\title{
A measure of Normalization Violence
}
\description{
Calculates 4 metrics to assess the degree of difference between normalized and raw betas.
 }
\usage{
qual(norm,raw)
}
\arguments{
  \item{norm}{
Matrix of normalized betas 
}
  \item{raw}{
Matrix of raw betas
}
}
%\details{}

\value{
Returns data.frame containing rmsd, sdd, sadd and srms for each sample (columns) in supplied matrices.
}

%\references{}

\examples{
library(wateRmelon)
data(melon)
d.melon <- dasen(melon)
raw.bet <- betas(melon)
norm.bet <- betas(d.melon)
qual(norm=norm.bet, raw=raw.bet)
}

\author{
Leo Schalkwyk
}

%\note{}

