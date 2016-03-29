\name{readEPIC}
\alias{readEPIC}
\title{
readEPIC
}
\description{
Reads Epic arrays from raw idats into \code{MethyLumiSet} objects from directory.
 }
\usage{
readEPIC(idatPath, barcodes=NULL, pdat=NULL,parallel=F,n=T,oob=F, ...) 
}
\arguments{
  \item{path}{
Path directory where idats are located. \code{readEPIC} looks in the specified path and
converts all .idats within path to relevant barcodes, which is then passed to a modified
version of \code{methylumIDAT} to parse all idats present in the specified directory. 
}
  \item{...}{
Additional arguments passed to \code{methylumIDAT}
}
}
\details{
Read a set of .idat files within a file directory and return a \code{MethylumiSet} object
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
}

\examples{
#Ficticious file pathway 
# path <- "Data/Experiment/Idatlocation"
# data <- readEPIC(path, pdat=NULL, oob=F, n=T)
}

