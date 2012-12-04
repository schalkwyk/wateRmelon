\name{dasen-minfi-methods}
\docType{methods}
\alias{naten,MethylSet-method}
\alias{naten,RGChannelSet-method}
\alias{betaqn,MethylSet-method}
\alias{betaqn,RGChannelSet-method}
\alias{nanet,MethylSet-method}
\alias{nanet,RGChannelSet-method}
\alias{nanes,MethylSet-method}
\alias{nanes,RGChannelSet-method}
\alias{danes,MethylSet-method}
\alias{danes,RGChannelSet-method}
\alias{danet,MethylSet-method}
\alias{danet,RGChannelSet-method}
\alias{daten1,MethylSet-method}
\alias{daten1,RGChannelSet-method}
\alias{daten2,MethylSet-method}
\alias{daten2,RGChannelSet-method}
\alias{nasen,MethylSet-method}
\alias{nasen,RGChannelSet-method}
\alias{dasen,MethylSet-method}
\alias{dasen,RGChannelSet-method}
\alias{danen,MethylSet-method}
\alias{danen,RGChannelSet-method}
\alias{tost,MethylSet-method}
\alias{tost,RGChannelSet-method}
\alias{fuks,MethylSet-method}
\alias{fuks,RGChannelSet-method}
\alias{swan,MethylSet-method}
\alias{swan,RGChannelSet-method}

\title{
Calculate normalized betas from Illumina 450K methylation arrays
}

\description{
Multiple ways of calculating the index of methylation (beta) from methylated and unmethylated probe intensities 
used in Pidsley et al 2012.  
}
\section{methods}{
\describe{
dasen ( mns, uns, onetwo, fudge = 100, ... )
nasen ( mns, uns, onetwo, fudge = 100 )
betaqn( bn ) 
naten ( mn,  un,          fudge = 100 )
naten ( mn,  un,          fudge = 100 )
nanet ( mn,  un,          fudge = 100 )
nanes ( mns, uns, onetwo, fudge = 100 )
danes ( mn,  un,  onetwo, fudge = 100, ... ) 
danet ( mn,  un,  onetwo, fudge = 100, ... )
danen ( mns, uns, onetwo, fudge = 100, ... )
daten1( mn,  un,  onetwo, fudge = 100, ... )
daten2( mn,  un,  onetwo, fudge = 100, ... )
tost  ( mn,  un,  da, pn )
fuks  ( data, anno)
swan  ( mn,  un,  qc )
}
}
\arguments{

  \item{mn, mns}{
matrix of methylated signal intensities, each column representing a sample 
(default method), or an object for which a method is available. Column names 
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

  \item{da, anno}{annotation data frame, such as x@featureData@data #methylumi package 
   }

  \item{qc}{
control probe intensities: list of 2 matrices, Cy3 and Cy5, with rownames, such
as produced by  intensitiesByChannel(QCdata(x)) #methylumi package 
   }

  \item{fudge}{
value added to total intensity to prevent denominators close to zero when
calculating betas
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
\bold{dasen}	same as nasen but type I and type II backgrounds are normalized first.  This is our recommended method

\bold{betaqn}	quantile normalizes betas 

\bold{naten}	quantile normalizes methylated and unmethylated intensities separately, then calculates betas

\bold{nanet}	quantile normalizes methylated and unmethylated intensities together, then calculates betas.  This should equalize dye bias.

\bold{nanes}	quantile normalizes methylated and unmethylated intensities separately, except for type II probes where methylated and unmethylated are normalized together. This should equalize dye bias without affecting type I probes which are not susceptible.

\bold{danes}	same as nanes, except typeI and type II background are equalised first.

\bold{danet}	same as nanet, except typeI and type II background are equalised first.

\bold{danen}	background equalisation only, no normalization

\bold{daten1}	same as naten,  except typeI and type II background are equalised first (smoothed only for methylated)

\bold{daten2}	same as naten,  except typeI and type II background are equalised first (smoothed for methylated an unmethylated)	

\bold{nasen}	same as naten but typeI and typeII intensities quantile normalized separately

\bold{tost}	method from Touleimat and Tost 2011

\bold{fuks}	method from Dedeurwaerder et al 2011.  Peak correction only, no normalization

\bold{swan}	method from Maksimovic et al 2012

}
\value{
a matrix of betas is returned by the MethySet and RGChannelSet methods because they do not have a defined slot for betas.

}
\references{
%% ~put references to the literature/web site here ~
[1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)

[2] Dedeurwaerder S, Defrance M, Calonne E, Sotiriou C, Fuks F: 
Evaluation of the Infinium Methylation 450K technology . Epigenetics 
2011, 3(6):771-784.

[3] Touleimat N, Tost J: Complete pipeline for Infinium R Human 
Methylation 450K BeadChip data processing using subset quantile 
normalization for accurate DNA methylation estimation. Epigenomics 
2012, 4:325-341)

[4] Maksimovic J, Gordon L, Oshlack A: SWAN: Subset quantile 
Within-Array Normalization for Illumina Infinium HumanMethylation450 
BeadChips. Genome biology 2012, 13(6):R44

}
\author{
%%  ~~who you are~~
Leonard.Schalkwyk@kcl.ac.uk
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
