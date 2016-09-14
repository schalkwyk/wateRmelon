\name{dasen}

\alias{dasen}
\alias{naten}
\alias{betaqn}
\alias{nanet}
\alias{nanes}
\alias{danes}
\alias{danet}
\alias{daten1}
\alias{daten2}
\alias{nasen}
\alias{danen}
\alias{tost}
\alias{fuks}
\alias{swan}

\title{
Calculate normalized betas from Illumina 450K methylation arrays
}

\description{
Multiple ways of calculating the index of methylation (beta) from 
methylated and unmethylated probe intensities used in Pidsley et al 
2012.  S4 methods exist where possible for MethyLumiSet, MethylSet, 
RGSet and exprmethy450 objects.
}

\usage{
dasen ( mns, uns, onetwo, fudge = 100, ret2=FALSE, ... )
nasen ( mns, uns, onetwo, ret2=FALSE, fudge = 100, ... )
betaqn( bn ) 
naten ( mn,  un, fudge = 100, ret2=FALSE, ... )
naten ( mn,  un, fudge = 100, ret2=FALSE, ... )
nanet ( mn,  un, fudge = 100, ret2=FALSE, ... )
nanes ( mns, uns, onetwo, fudge = 100, ret2=FALSE, ... )
danes ( mn,  un,  onetwo, fudge = 100, ret2=FALSE, ... ) 
danet ( mn,  un,  onetwo, fudge = 100, ret2=FALSE, ... )
danen ( mns, uns, onetwo, fudge = 100, ret2=FALSE, ... )
daten1( mn,  un,  onetwo, fudge = 100, ret2=FALSE, ... )
daten2( mn,  un,  onetwo, fudge = 100, ret2=FALSE, ... )
tost  ( mn,  un,  da, pn )
fuks  ( data, anno)
swan  ( mn,  un,  qc, da=NULL, return.MethylSet=FALSE )
}

\arguments{

  \item{mn, mns}{
matrix of methylated signal intensities, each column representing a sample 
(generic) or a MethyLumiSet, RGSet, or MethylSet object. Column names 
are used to get Sentrix row and column by default, see '...'.
   }

  \item{un, uns}{
matrix of unmethylated signal intensities, each column representing a sample
(default method) or NULL when mn is an object containing methylated and 
unmethylated values 
   }

  \item{bn, data}{
matrix of precalculated betas, each column representing a sample 
   }

  \item{onetwo}{
character vector or factor of length nrow(mn) indicating assay type 'I'
or 'II' 
   }

  \item{pn}{
matrix of detection p-values, each column representing a sample
   }

  \item{da, anno}{annotation data frame, such as x@featureData@data #methylumi package.  If NULL, the swan method requires the \code{IlluminaHumanMethylation450kmanifest} package. 
   }

  \item{qc}{
control probe intensities: list of 2 matrices, Cy3 and Cy5, with rownames, such
as produced by  intensitiesByChannel(QCdata(x)) #methylumi package 
   }

  \item{fudge}{
value added to total intensity to prevent denominators close to zero when
calculating betas
   }

  \item{return.MethylSet}{
if TRUE, returns a MethylSet object instead of a naked matrix of betas.
   }

  \item{ret2}{
if TRUE, returns a list of intensities and betas instead of a naked 
matrix of betas.
   }

  \item{...}{
additional argument roco for dfsfit giving Sentrix rows and columns. 
This allows a background gradient model to be fit.  This is split from 
data column names by default.  roco=NULL disables model fitting
(and speeds up processing), otherwise roco can be supplied as a 
character vector of strings like 'R01C01' (only 3rd and 6th characters used).
   }
}

\details{
\bold{dasen}	same as nasen but type I and type II backgrounds are equalized first.  This is our recommended method

\bold{betaqn}	quantile normalizes betas 

\bold{naten}	quantile normalizes methylated and unmethylated intensities separately, then calculates betas

\bold{nanet}	quantile normalizes methylated and unmethylated intensities together, then calculates betas.  This should equalize dye bias

\bold{nanes}	quantile normalizes methylated and unmethylated intensities separately, except for type II probes where methylated and unmethylated are normalized together. This should equalize dye bias without affecting type I probes which are not susceptible

\bold{danes}	same as nanes, except type I and type II background are equalized first

\bold{danet}	same as nanet, except type I and type II background are equalized first

\bold{danen}	background equalization only, no normalization

\bold{daten1}	same as naten,  except type I and type II background are equalized first (smoothed only for methylated)

\bold{daten2}	same as naten,  except type I and type II background are equalized first (smoothed for methylated an unmethylated)	

\bold{nasen}	same as naten but type I and typeII intensities quantile normalized separately

\bold{tost}	method from Touleimat and Tost 2011

\bold{fuks}	method from Dedeurwaerder et al 2011.  Peak correction only, no normalization

\bold{swan}	method from Maksimovic et al 2012

}
\value{
a matrix (default method) or object of the same shape and order as the 
first argument containing betas. 

}
\references{
[1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)

[2] Dedeurwaerder S, Defrance M, Calonne E, Sotiriou C, Fuks F: 
Evaluation of the Infinium Methylation 450K technology . Epigenetics 
2011, 3(6):771-784.

[3] Touleimat N, Tost J: Complete pipeline for Infinium R Human 
Methylation 450K BeadChip data processing using subset quantile 
normalization for accurate DNA methylation estimation. Epigenomics 
2012, 4:325-341.

[4] Maksimovic J, Gordon L, Oshlack A: SWAN: Subset quantile 
Within-Array Normalization for Illumina Infinium HumanMethylation450 
BeadChips. Genome biology 2012, 13(6):R44

}
\author{
Leonard.Schalkwyk@kcl.ac.uk
}
%\note{


\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{pfilter}}, \code{\link{as.methylumi}}
}

\examples{

#MethyLumiSet method
data(melon)
melon.dasen <- dasen(melon)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
