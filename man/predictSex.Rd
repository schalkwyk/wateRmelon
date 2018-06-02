\name{predictSex}
\alias{predictSex}

\title{
Sex Prediction of samples assayed on the Illumina Infinium BeadChip platform
}

\description{
Predicts the sex of samples
}

\usage{
predictSex(x, x.probes, pc = 2, plot=FALSE, irlba=TRUE, center=FALSE, scale. = FALSE)
}

\arguments{
\item{x}{A matrix containing values equivalent to betas.
}
\item{x.probes}{vector containing logical or numeric index of probes on the X-chromosome
}
\item{pc}{The principal component to guess sex across (usually the 2nd one in most cases)
}
\item{plot}{Logical, indicated whether or not to plot the prediction}
\item{irlba}{Logical, indicate whether or not to use the faster method to generate principal components}
\item{center}{Logical, indicate whether or not to center data around 0}
\item{scale.}{Logical, indicate whether or not to scale data prior to prcomp}
}

\value{
Character Vector of predicted sex for each sample, samples that were not guess are assigned 'Unknown'
}

\details{
Function attempts to predictSex based on the colMeans of the X chromsome and princpal components. Assumes that there are two sexes to predict and as a result will perform poorly on data-sets consisting of a single data-set.
}

\author{
Tyler J. Gorrie-Stone
Who to contact: <tgorri@essex.ac.uk>
}

%\examples{

