\name{pwod-methods}
\docType{methods}
\alias{pwod,MethylSet-method}
\alias{pwod,MethyLumiSet-method}
\alias{pwod,RGChannelSet-method}
\title{
Methods for Function \code{pwod} in Package \pkg{wateRmelon}
}
\description{
Methods for function \code{pwod}, please see \code{pwod} for details
of how function performs.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "MethyLumiSet")}}{
all of the methods simply extract betas from the data object (which can be
a \code{MethylSet}, \code{MethyLumiSet}, or 
\code{RGChannelSet}) and calculates the outliers.
}
}}
\keyword{methods}