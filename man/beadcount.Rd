\name{beadcount}
\alias{beadcount}
\title{
Creates matrix of beacounts from minfi data.
}
\description{
Creates matrix of beacounts from data read in using the minfi package.  NAs represent probes with beadcount <3.  An Extended RG Channel Set is required for this function to work.
}
\usage{
beadcount(x)
}

\arguments{
  \item{x}{
450K methylation data read in using minfi to create an Extended RG Channel Set
}
}
\value{
A matrix of bead counts with bead counts <3 represented by NA for use in the pfilter function for quality control
}
\references{
[1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
}
\author{
Ruth.Pidsley@kcl.ac.uk
}
\note{
The beadcount function is internal to the pfilter function
}
