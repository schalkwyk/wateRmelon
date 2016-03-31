\name{got}
\alias{got}
\alias{fot}
%\alias{pop}
\title{
Internal functions for Illumina i450 normalization functions
}
\description{
got and fot find the annotation column differentiating type I and type II 
assays in MethylSet (got) or MethyLumiSet (fot) objects. pop extracts columns
from IlluminaHumanMethylation450k.db
}
\usage{
got(obj)
fot(x)
%pop(fd,rn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  a MethyLumiSet }
  \item{obj}{  a MethylSet  }
%  \item{fd}{ a character vector of the desired annotation columns    }
%  \item{rn}{ a character vector of the desired features   }

}
\details{
\code{got} returns a character vector of 'I' and 'II', \code{fot} returns the index of the relevant column.
\code{pop} returns a data frame
}
\references{
Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)}
\author{
Leonard.Schalkwyk@kcl.ac.uk
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
