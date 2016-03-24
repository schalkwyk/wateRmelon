\name{bscon}
\alias{bscon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate bisulphite conversion
%%  ~~function to do ... ~~
}
\description{
Uses control data from Infinium HumanMethylation450 BeadChip Kit to calculate bisulfite conversion for each array
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bscon(x, ...)  # S4 methods exist for RGset and MethyLumiset objects
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ IDAT or report files containing 450k data 
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{ current methods have no optional arguments
%%     ~~Describe \code{\dots} here~~
}
}
\details{ This function uses the green and red channels reading of the type I and type II bisilfite conversion data to return the median bisulfite conversion percentage value for each array.

For the type I chemistry the beta values are calculated by dividing the first three probes of the green channel (C1, C2, C3) and the second three probes of the red channel (C4, C5, C6) by the sum of these probes and the unconverted probes of the green (U1, U2, U3) and the red (U4, U5, U6) channel.

The beta values from type II chemistry are calculated by dividing the methylated (red) channels by the sum of methylated (red) and unmethylated (green) channels.


%%  ~~ If necessary, more details than the description above ~~
}
\value{A vector of percentage values referring to the bisulfite conversion levels of each array.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Louis El Khoury (louis.el-khoury@essex.ac.uk), Eilis Hannon, Leonard Schalkwyk (lschal@essex.ac.uk)
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(wateRmelon)
data(melon)
bs <- bscon(melon)
bs
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QC data }
\keyword{ Bisulphite Conversion Rate}% __ONLY ONE__ keyword per line
