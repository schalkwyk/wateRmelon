\name{estimateCellCounts}
\alias{estimateCellCounts.wmln}
\alias{estimateCellCounts.wateRmelon}
\alias{estimateCellCounts.wateRmelon,MethylSet-method}
\alias{estimateCellCounts.wateRmelon,MethyLumiSet-method}
\alias{estimateCellCounts.wateRmelon,RGChannelSet-method}
\title{
Cell Proportion Estimation using wateRmelon
}
\description{
Estimates relative proportion of pure cell types within a sample, mostly identical to \code{\link[minfi]{estimateCellCounts}}. 
References for both 450k and EPIC array are available. However 450k reference can be used on EPIC data by specifying the reference platform. Additionally a measure of error is calculated as a means of quality control.
}
\usage{
    estimateCellCounts.wmln(
        object,
        referencePlatform = c("IlluminaHumanMethylation450k",
            "IlluminaHumanMethylationEPIC",
            "IlluminaHumanMethylation27k"),
        mn = NULL,
        un = NULL,
        bn = NULL,
        perc = 1,
        compositeCellType = "Blood",
        probeSelect = "auto",
        cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
        returnAll = FALSE,
        meanPlot = FALSE,
        verbose=TRUE,
        ...)
}
\arguments{
\item{object}{An object of class methylumiset, which contains (un)normalised methylated and unmethylated intensities}
\item{mn}{
if NULL will call methylated(object), otherwise can be given matrix of identical dimensions to object. 
}
\item{un}{if NULL will call unmethylated(object), otherwise can be given matrix of identical dimensions to object.}

\item{bn}{if NULL will call betas(object), otherwise can be given matrix of identical dimensions to object.}
\item{perc}{Percentage of query-samples to use to normalise reference dataset. This should be 1 unless using a very large data-set then lowering this will allow for an increase in performance}
\item{compositeCellType}{Which composite cell type is being deconvoluted. Should be either "Blood", "CordBlood", or "DLPFC"}
  \item{probeSelect}{How should probes be selected to distinguish cell types? Options include
    "both", which selects an equal number (50) of probes (with F-stat p-value < 1E-8) with the
    greatest magnitude of effect from the hyper- and hypo-methylated sides, and "any", which
    selects the 100 probes (with F-stat p-value < 1E-8) with the greatest magnitude of difference
    regardless of direction of effect. Default input "auto" will use "any" for cord blood and
    "both" otherwise, in line with previous versions of this function and/or our recommendations.
    Please see the references for more details.}
  \item{cellTypes}{Which cell types, from the reference object, should be
    we use for the deconvolution? See details. }
  \item{referencePlatform}{The platform for the reference dataset; if
    the input \code{rgSet} belongs to another platform, it will be
    converted using \code{\link[minfi]{convertArray}}.}
  \item{returnAll}{Should the composition table and the normalized user
    supplied data be return?}
  \item{verbose}{Should the function be verbose?}
  \item{meanPlot}{
    Whether to plots the average DNA methylation across the cell-type discrimating probes within
    the mixed and sorted samples.
}
\item{...}{Other arguments, i.e arguments passed to plots}
}
\details{
See \code{\link[minfi]{estimateCellCounts}} for more information regarding the exact details. estimateCellCounts.wmln differs slightly, as it will impose the quantiles of type I and II probes onto the reference Dataset rather than normalising the two together. This is 1) More memory efficient and 2) Faster - due to not having to normalise out a very small effect the other 60 samples from the reference set will have on the remaining quantiles.

Optionally, a proportion of samples can be used to derive quantiles when there are more than 1000 samples in a dataset, this will further increase performance of the code at a cost of precision. If data is pre-normalised a minimum of two samples are required.
}
