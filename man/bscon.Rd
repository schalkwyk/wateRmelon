\name{bscon}
\alias{bscon}
\alias{bscon,MethyLumiSet-method}
\alias{bscon,RGChannelSet-method}
\title{ Calculate bisulphite conversion
}
\description{
Uses control data from Infinium HumanMethylation450 BeadChip to calculate bisulfite conversion for each array
}
\usage{
bscon(x, ...)  # S4 methods exist for RGChannelSet and MethyLumiset objects
}
\arguments{
  \item{x}{ IDAT or report files containing 450k data 
}
  \item{\dots}{ current methods have no optional arguments
}
}
\details{ This function uses the green and red channels reading of the type I and type II bisilfite conversion data to return the median bisulfite conversion percentage value for each array.

For the type I chemistry the beta values are calculated by dividing the first three probes of the green channel (C1, C2, C3) and the second three probes of the red channel (C4, C5, C6) by the sum of these probes and the unconverted probes of the green (U1, U2, U3) and the red (U4, U5, U6) channel.

The beta values from type II chemistry are calculated by dividing the methylated (red) channels by the sum of methylated (red) and unmethylated (green) channels.

}

\value{A vector of percentage values referring to the bisulfite conversion levels of each array.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{}
\author{Louis El Khoury (louis.el-khoury@essex.ac.uk), Eilis Hannon, Leonard Schalkwyk (lschal@essex.ac.uk)
}

\note{
Updates to HumanMethylationEPIC manifest has seen the removal of control probes C6 and U6. 
This does not appear to grossly affect how function performs however we are 
considering alternative approaches to account for this.
}


\examples{
library(wateRmelon)
data(melon)
bs <- bscon(melon)
bs
}
\keyword{ QC data }
\keyword{ Bisulphite Conversion Rate}% __ONLY ONE__ keyword per line
