\name{adaptRefQuantiles}
\alias{adaptRefQuantiles}
\alias{coRankedMatrices}
\alias{concatenateMatrices}
\alias{dataDetectPval2NA}
\alias{detectionPval.filter}
\alias{filterXY}
\alias{findAnnotationProbes}
\alias{getMethylumiBeta}
\alias{getQuantiles}
\alias{getSamples}
\alias{loadMethylumi2}
\alias{lumiMethyR2}
\alias{nbBeadsFilter}
\alias{normalize.quantiles2}
\alias{normalizeIlluminaMethylation}
\alias{pipelineIlluminaMethylation.batch}
\alias{preprocessIlluminaMethylation}
\alias{referenceQuantiles}
\alias{robustQuantileNorm_Illumina450K}
\alias{robustQuantileNorm_Illumina450K.probeCategories}
\alias{uniqueAnnotationCategory}
\alias{bgIntensitySwan.methylumi}

\title{ Functions from 450-pipeline  (Touleimat & Tost)
}
\description{
These functions are part of the 450K pipeline (Touleimat and Tost, Epigenomics 2012 4:325).  For freestanding use of the normalization function, a wrapper is 
provided, see \code{\link{tost}}
}
%\usage{

%tost(mn,  un,  da, pn )

%}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%  \item{mn}{
%matrix of methylated signal intensities, each column representing a sample 
%(default method), or an object for which a method is available. Column names 
%are used to get Sentrix row and column by default, see '...'.
%   }
%
%  \item{un}{
%matrix of unmethylated signal intensities, each column representing a sample
%(default method) or NULL when mn is an object containing methylated and 
%unmethylated values 
%   }
%
%  \item{da}{annotation data frame, such as x@featureData@data #methylumi package 
%   }
%
%
%}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%see \code{\link{tost}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%
%}
\references{
Touleimat N, Tost J: Complete pipeline for Infinium R Human 
Methylation 450K BeadChip data processing using subset quantile 
normalization for accurate DNA methylation estimation. Epigenomics 
2012, 4:325-341

Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
}
\author{
Nizar Touleimat, wrapper by Leonard.Schalkwyl@kcl.ac.uk
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}


