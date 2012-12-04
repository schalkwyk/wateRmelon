\name{beadc}
\alias{beadc}
\title{
Calculates the number of samples with bead count <3 for each probe in matrix of bead count values
}
\description{
Calculates the number of samples with bead count <3 for each probe in matrix of bead count values.
 }
\usage{
beadc(x)
}
\arguments{
  \item{x}{
matrix of bead count values returned by the beadcount function
}
}
\value{
Vector of number of samples with bead count <3 for each probe
}
\references{
[1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)}
\author{
ruth.pidsley@kcl.ac.uk
}
\note{
The beadc function is internal to the pfilter function
}

