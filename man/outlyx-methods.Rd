\name{outlyx-methods}
\docType{methods}
\alias{outlyx,MethylSet-method}
\alias{outlyx,RGChannelSet-method}
\alias{outlyx,MethyLumiSet-method}
\title{
Methods for Function \code{outlyx} in Package \pkg{wateRmelon}
}
\description{
Methods for function \code{outlyx}, please see \code{outlyx} for details
of how function performs.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "MethyLumiSet")}}{
all of the methods simply extract betas from the data object (which can be
a \code{MethylSet}, \code{MethyLumiSet}, or 
\code{RGChannelSet}) and calculates the outliers.
}
}}
\keyword{methods}