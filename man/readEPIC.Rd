\name{readEPIC}
\alias{readEPIC}
\title{
readEPIC
}
\description{
Reads Epic arrays from raw idats into \code{MethyLumiSet} objects from directory.
 }
\usage{
#Ficticious file pathway 
path <- "Data/Experiment/Idatlocation"
data <- readEPIC(path, ...)
}
\arguments{
  \item{path}{
Path directory where idats are located. \code{readEPIC} uses string manipulation
to generate a vector of barcodes which is then passed to a modified
version of \code{methylumIDAT} to read in all idats present in the directory. 
}
  \item{...}{
See \code{methylumIDAT}
}
}
\details{
Read a set of .idat files within a file directory and return a \code{MethylumiSet}
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

If barcodes are the preffered input, use \code{methylumIDAT} as normal. 
Specifiying idatPath.
}

