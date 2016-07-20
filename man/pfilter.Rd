\name{pfilter}
\alias{pfilter-methods}
\alias{pfilter,MethyLumiSet-method}
\alias{pfilter,RGChannelSetExtended-method}
\alias{pfilter}
\title{
  Basic data filtering for Illumina 450 methylation data 
}
\description{
The pfilter function filters data sets based on bead count and 
detection p-values.  The user can set their own thresholds or 
use the default pfilter settings.  pfilter will take data matrices 
of beta values, signal intensities and annotation data, but will 
also take methylumi (MethyLumiSet) or minfi (RGChannelSetExtended) 
objects.  However it has come to our attention that data read in 
using the various packages and input methods will give subtly variable
data output as they calculate detection p-value and beta values differently, 
and do/don?t give information about beadcount.  The pfilter function does not correct 
for this, but simply uses the detection p-value and bead count 
provided by each package.
}
\usage{
pfilter(mn, un, bn, da, pn, bc, perCount=NULL, pnthresh = NULL, perc = NULL, 
pthresh = NULL,logical.return=FALSE)
}


\section{Methods}{
   \describe{
      \item{\code{signature(mn = "MethyLumiSet")}}{
This is used for performing the pfilter method on MethyLumiSet objects produced by 
methylumiR.
      }
      \item{\code{signature(mn = "RGChannelSetExtended")}}{
This is used for performing the pfilter method on RGChannelSetExtended 
objects produced by minfi.
      }
   }
}




\arguments{
  \item{mn}{
matrix of methylated signal intensities, each column representing a sample 
(default method), or an object for which a method is available e.g 
MethyLumiSet or RGChannelSetExtended.    N.B. Bead count filtering will 
not work unless data read in as an  RGGhannelSetExtended rather than an 
RGChannelSet. 
}
  \item{un}{
matrix of unmethylated signal intensities, each column representing a sample
(default method) or NULL when mn is a MethyLumiSet or RGChannelSetExtended 
object 
}
  \item{bn}{
matrix of precalculated betas, each column representing a sample, or NULL when 
mn is a MethyLumiSet or RGChannelSetExtended object 
}
  \item{da}{
annotation data frame, such as x@featureData@data #methylumi package, or 
NULL when mn is a MethyLumiSet or RGChannelSetExtended object 

}
  \item{pn}{
matrix of detection p-values, each column representing a sample, a 
MethyLumiSet or RGChannelSetExtended object 
}
\item{bc}{
matrix of arbitrary values, each column representing a sample and eeach row 
representing a probe, in which "NA" represents beadcount <3, or NULL when mn 
is a MethyLumiSet or RGChannelSetExtended object 
}
  \item{perCount}{
remove sites having this percentage of samples with a beadcount <3, default set to 5
}

  \item{pnthresh}{
cutoff for detection p-value, default set to 0.05
}
  \item{perc}{
remove samples having this percentage of sites with a detection p-value 
greater than pnthresh, default set to 1
}
  \item{pthresh}{
remove sites having this percentage of samples with a detection p-value 
greater than pnthresh, default set to 1
}

  \item{logical.return}{
If it is TRUE, FALSE or TRUE is returned to indicate success
}
}
\value{
a filtered MethyLumiSet or

a list of the filtered matrices:

mn : methylated intensities

un : unmethylated intensities

bn : betas

da : feature data

or

a filtered MethylSet object.
}
\references{
[1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
}
\author{
Ruth.Pidsley@kcl.ac.uk
}

\note{
Adjusted \code{RGChannelSetExtended} methods - 12/10/2015
Now outputs a MethylSet object using preprocessRaw from minfi.
}

\examples{
# MethyLumiSet method
data(melon)
melon.pf <- pfilter(melon)
}
