\name{dasen-methods}

\alias{naten,MethyLumiSet-method}
\alias{betaqn,MethyLumiSet-method}
\alias{nanet,MethyLumiSet-method}
\alias{nanes,MethyLumiSet-method}
\alias{danes,MethyLumiSet-method}
\alias{danet,MethyLumiSet-method}
\alias{daten1,MethyLumiSet-method}
\alias{daten2,MethyLumiSet-method}
\alias{nasen,MethyLumiSet-method}
\alias{dasen,MethyLumiSet-method}
\alias{danen,MethyLumiSet-method}
\alias{tost,MethyLumiSet-method}
\alias{fuks,MethyLumiSet-method}
\alias{swan,MethyLumiSet-method}

\title{
Calculate normalized betas from MethyLumiSets of Illumina 450K methylation arrays
}
\description{
Multiple ways of calculating the index of methylation (beta) from methylated and unmethylated probe intensities 
used in Pidsley et al 2012.  
}

\section{methods}{
\describe{
dasen ( mns, fudge = 100, ... )
nasen ( mns, fudge = 100 )
betaqn( bn ) 
naten ( mn, fudge = 100 )
naten ( mn, fudge = 100 )
nanet ( mn, fudge = 100 )
nanes ( mns,fudge = 100 )
danes ( mn, fudge = 100, ... ) 
danet ( mn, fudge = 100, ... )
danen ( mns,fudge = 100, ... )
daten1( mn, fudge = 100, ... )
daten2( mn, fudge = 100, ... )
tost  ( mn )
fuks  ( data)
swan  ( mn )
}
}


\arguments{

  \item{mn, mns, data, bn}{
A MethyLumiSet object. Sample names names 
are used to get Sentrix row and column by default, see '...'.
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
a matrix (default method) or object of the same shape and order as the 
first argument containing betas. 

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
2012, 4:325-341

[4] Maksimovic J, Gordon L, Oshlack A: SWAN: Subset quantile 
Within-Array Normalization for Illumina Infinium HumanMethylation450 
BeadChips. Genome biology 2012, 13(6):R44

}
\author{
Leonard.Schalkwyk@kcl.ac.uk
}


