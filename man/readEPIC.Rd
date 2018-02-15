\name{readEPIC}
\alias{readEPIC}
\alias{epic.controls}
\title{
readEPIC
}
\description{
Reads Epic arrays from raw idats into \code{MethyLumiSet} objects from directory.
 }
\usage{
readEPIC(idatPath, barcodes=NULL, pdat=NULL,parallel=F,n=T,oob=F,force=F, ...) 
}
\arguments{
  \item{idatPath}{
Path directory where .idat files are located. \code{readEPIC} looks in the specified path and
converts all .idats within path to relevant barcodes, which is then passed to a modified
version of \code{methylumIDAT} to parse all idats present in the specified directory.
}
  \item{barcodes}{
If NULL, function will search supplied argument in "idatPath" for all idats within directory. 
If given a vector of barcodes, readEPIC will search for those specific barcodes within the idatPath supplied. 
}
  \item{parallel}{If TRUE, an attempt will be made to process using multiple
          cores on a multicore machine.
}
  \item{pdat}{A data.frame describing the samples.  A special column named
          "barcodes" can be used to specify the barcodes to be read when using \code{methylumIDATepic}. See \code{methylumIDAT for usage}
}
  \item{n}{If TRUE, beadcounts from .idats will be included in final object}
  \item{oob}{If TRUE, out-of-band (OOB) or opposite-channel signals will be kept}
  \item{force}{If TRUE, will combine EPIC IDATs read with differing dmaps}
  \item{...}{
Additional arguments passed to \code{methylumIDAT}
}
}
\details{
Read a set of .idat files within a file directory and return a \code{MethylumiSet} object. 
}

\value{
A \code{MethyLumiSet} object.
}
\references{
methylumi
}

\author{
Tyler Gorrie-Stone - tgorri@essex.ac.uk
}

\note{
Contains heavily modified version of \code{methylumIDAT} and other accessory 
functions used to construct a \code{MethylumiSet} object, specifically tailored
for EPIC arrays. \code{readEPIC} can also handle 450k and 27k arrays as 
\code{methylumIDAT} functionality for these platforms remains unchanged. 

Alternatively it is possible to invoke \code{methylumIDATepic} to use the modified
version \code{methylumIDAT}, which has similar usage.

EPIC manifest has since been updated to B4 version, which has notably fewer probes than previous manifests. It is entirely possible that we will migrate to the manifest packages available on BioConductor and allow for versioning control.
}

\examples{
#Ficticious file pathway 
# path <- "Data/Experiment/Idatlocation"
# data <- readEPIC(path, barcodes = NULL oob=F, n=T)
}

