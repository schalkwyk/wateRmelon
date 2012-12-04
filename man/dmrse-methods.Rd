\name{dmrse-methods}
\docType{methods}
\alias{dmrse-methods}
\alias{dmrse,exprmethy450-method}
\alias{dmrse,MethylSet-method}
\alias{dmrse,MethyLumiSet-method}
\alias{dmrse,RGChannelSet-method}
\alias{dmrse_col-methods}
\alias{dmrse_col,exprmethy450-method}
\alias{dmrse_col,MethylSet-method}
\alias{dmrse_col,MethyLumiSet-method}
\alias{dmrse_col,RGChannelSet-method}
\alias{dmrse_row-methods}
\alias{dmrse_row,exprmethy450-method}
\alias{dmrse_row,MethylSet-method}
\alias{dmrse_row,MethyLumiSet-method}
\alias{dmrse_row,RGChannelSet-method}
\title{ Methods for Function \code{dmrse} in Package \pkg{wateRmelon} }
\description{
Methods for function \code{dmrse}, \code{dmrse_row} and 
\code{dmrse_col} in package \pkg{wateRmelon}.  Please see 
\code{\link{dmrse}} for details of the calculation of this 
standard-error performance metric.
}
\section{Methods}{
\describe{


\item{\code{signature(betas = "exprmethy450")}}{
 all of the methods simply extract betas from the data object (which can be
a \code{exprmethy450}, \code{MethylSet}, \code{MethyLumiSe}, or 
\code{RGChannelSet}) and calculate the metric.
}

}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
