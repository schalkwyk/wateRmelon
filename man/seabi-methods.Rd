\name{seabi-methods}
\docType{methods}
\alias{seabi-methods}
\alias{seabi,exprmethy450-method}
\alias{seabi,MethylSet-method}
\alias{seabi,MethyLumiSet-method}
\alias{seabi,RGChannelSet-method}
\title{Methods for Function \code{seabi} in Package \pkg{wateRmelon}}
\description{
Methods for function \code{seabi} in package \pkg{wateRmelon}. 
Please see \code{\link{seabi}} for details of the calculation of 
this ROC AUC performance metric. 
}
\section{Methods}{
\describe{
\item{\code{signature(betas = "exprmethy450")}}{
 all of the methods simply extract betas from the data object 
(which can be a \code{exprmethy450}, \code{MethylSet}, 
\code{MethyLumiSe}, or \code{RGChannelSet}) and calculate 
the metric. All the methods also require a factor differentiating 
male from female samples.
}

}

}
\keyword{methods}
