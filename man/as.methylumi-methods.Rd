\name{as.methylumi-methods}
%\docType{methods}
\alias{as.methylumi}
\alias{as.methylumi-methods}
\alias{as.methylumi,ANY-method}
\alias{as.methylumi,MethylSet-method}
\alias{as.methylumi,MethyLumiSet-method}
\alias{getColumns}
\title{  Methods for Function \code{as.methylumi}  }
\description{
Returns a MethyLumiSet object populaed with the data provided.  There are
MethyLumiSet and MethylSet methods.  In the default method, the data
is all optional.  Please note that for the results to be
sane, mn, un, bn, and pv have to be in the same sample and
feature order and the same size.  The function does not currently
do any checks!
}

\usage{

# default method
as.methylumi (mn = NULL, un = NULL, bn = NULL, pv = NULL, qc = NULL, da = NULL,...)

}


\section{Methods}{
   \describe{
      \item{\code{signature(mn = "MethylSet")}}{
Coerces a MethylSet to a MethyLumiSet, and provides it with a set of
featureData, which by default is just the chromosome and DESIGN (ie typeI
or type II assay).  Other data can be included using the \code{fd}
argument, available data is listed by the function \code{getColumns()}
      }

      \item{\code{signature(mn = "MethyLumiSet")}}{
This is mainly useful for adding featureData as described under MethylSet above.
MethyLumiSet objects produced by methylumiR have the full annotation, those
from methylumIDAT do not, and functions such as \code{\link{swan}} require it

      }
      \item{\code{signature(mn = "ANY")}}{
\code{as.methylumi (mn = NULL, un = NULL, bn = NULL, pv = NULL, qc = NULL,
    da = NULL, fd = c("CHR", "DESIGN")) }
      }
   }
}


\arguments{

  \item{mn}{
matrix of methylated signal intensities, each column representing a sample
(generic) or a MethyLumiSet, RGSet, or MethylSet object. Column names
are used to get Sentrix row and column by default, see '...'.
   }

  \item{un}{
matrix of unmethylated signal intensities, each column representing a sample
(default method) or NULL when mn is an object containing methylated and
unmethylated values
   }

  \item{bn}{
matrix of precalculated betas, each column representing a sample
   }

  \item{pv}{
matrix of detection p-values, each column representing a sample
   }

  \item{da}{annotation data frame, such as x@featureData@data #methylumi package.  If NULL (the default), the  \code{IlluminaHumanMethylation450kmanifest} package is used.  See the fd argument
   }

  \item{qc}{
control probe intensities: list of 2 matrices, Cy3 and Cy5, with
rownames, such as produced by  intensitiesByChannel(QCdata(x))
(methylumi package)
   }
%  \item{fd}{  vector of items of featureData, which by default is
%just the chromosome and DESIGN (ie typeI or type II assay).  Other
%data can be included using the \code{fd} argument, available data
%is listed by the function \code{getColumns()}
%   }
%  \item{ad}{ optional assayData }
\item{...}{Other arguments such as a featureData object or optional assayData}
}


\keyword{methods}
