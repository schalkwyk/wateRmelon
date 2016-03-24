\name{combo}
\alias{combo}
\title{ Combine MethyLumiSet objects
}
\description{
This is a wrapper for combining different MethyLumiSet objects.
}
\usage{
combo(...)
}
\arguments{
  \item{\dots}{
Eventually, any number of MethyLumiSet objects. Currently only guaranteed for 
2 objects. 
}
}
\details{
This is a wrapper for \code{methylumi::combine}, which works around a name clash with
a different combine function from the \code{gdata} package, and also a bug in 
\code{methylumi::combine}. 
}
\value{
a \code{MethyLumiSet}.  The \code{assayData}, \code{QCdata}, \code{experimentData}, \code{protocolData}
and \code{phenoData} are joined on \code{sampleName} .  \code{featureData}  and {annotation} are taken from
the object given in the first argument

}
\references{
[1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)}

\author{
Leo Schalkwyk <leonard.schalkwyk@kcl.ac.uk>

}
\note{
the function uses \code{sampleNames} and gets rid of duplicates. 
Numeric sampleNames cause problems (and are a Bad Idea anyway).  They
should be turned into names with \code{make.names()} first.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{as.methylumi}}
}
\examples{

library(wateRmelon)
data(melon)
## pretend we have two different data sets
melon
pelon <- melon
sampleNames(pelon) <- gsub('^6', 7, sampleNames(pelon))
combo(melon, pelon)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
