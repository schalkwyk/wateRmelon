\name{internal}
\alias{DataToNChannelSet2}
\alias{IDATsToMatrices2}
\alias{IDATtoMatrix2}
\alias{NChannelSetToMethyLumiSet2}
\alias{anSNP}
\alias{anti.trafo}
\alias{aoget}
\alias{bfp}
\alias{bscon_methy
\alias{bscon_minfi}
\alias{columnMatrix}
\alias{designIItoMandU2}
\alias{designItoMandU2}
\alias{dfort}
\alias{extractAssayDataFromList2}
\alias{genall}
\alias{genme}
\alias{genus}
\alias{getControlProbes2}
\alias{getMethylationBeadMappers2}
\alias{goodSNP}
\alias{iqrFun}
\alias{mergeProbeDesigns2}
\alias{methylumIDATepic}
\alias{mvFun}
\alias{oxyscale}
\alias{pcouted}
\alias{trafo}
\title{
Internal functions for readEPIC and other wateRmelon functions introduced in v 1.13.1
}
\description{
few if any functions of interest to users
}
\usage{
DataToNChannelSet2(mats, chans = c(Cy3 = "GRN", Cy5 = "RED"), parallel = F, protocol.data = F, IDAT = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mats}{
%%     ~~Describe \code{mats} here~~
}
  \item{chans}{
%%     ~~Describe \code{chans} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{protocol.data}{
%%     ~~Describe \code{protocol.data} here~~
}
  \item{IDAT}{
%%     ~~Describe \code{IDAT} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mats, chans = c(Cy3 = "GRN", Cy5 = "RED"), parallel = F, 
    protocol.data = F, IDAT = TRUE) 
{
    epic = hm27 = hm450 = 0
    qw <- unlist(lapply(mats, function(x) attr(x, "ChipType")))
    epic = sum(grepl("BeadChip 8x5", qw))
    message(paste(epic, "Epic samples found"))
    hm450 = sum(grepl("BeadChip 12x8", qw))
    message(paste(hm450, "HumanMethylation450 samples found"))
    hm27 = sum(grepl("BeadChip 12x1", qw))
    message(paste(hm27, "HumanMethylation27 samples found"))
    if (hm27 > 0 && hm450 > 0 | hm27 > 0 && epic > 0 | hm450 > 
        0 && epic > 0) {
        stop("Cannot process both platforms simultaneously; please run separately.")
    }
    stopifnot(is(mats, "list"))
    assayNames = paste0(names(chans), c(".Mean"))
    assayNames = c(assayNames, paste0(names(chans), c(".NBeads")))
    names(assayNames) = assayNames
    fnames <- rownames(mats[[1]])
    extract <- function(assay) extractAssayDataFromList2(assay, 
        mats, fnames)
    assays = lapply(assayNames, extract)
    nb <- pmin(assays[["Cy3.NBeads"]], assays[["Cy5.NBeads"]])
    obj = new("NChannelSet", assayData = assayDataNew(R = assays[["Cy5.Mean"]], 
        G = assays[["Cy3.Mean"]], N = nb))
    featureNames(obj) = fnames
    if (IDAT) {
        message("Attempting to extract protocolData() from list...")
        ChipType = attr(mats[[1]], "ChipType")
        RunInfo = lapply(mats, function(d) attr(d, "RunInfo"))
        if (protocol.data) {
            scanDates = data.frame(DecodeDate = rep(NA, length(mats)), 
                ScanDate = rep(NA, length(mats)))
            rownames(scanDates) = names(mats)
            for (i in seq_along(mats)) {
                cat("decoding protocolData for", names(mats)[i], 
                  "...\n")
                if (nrow(RunInfo[[i]]) >= 2) {
                  scanDates$DecodeDate[i] = RunInfo[[i]][1, 1]
                  scanDates$ScanDate[i] = RunInfo[[i]][2, 1]
                }
            }
            protocoldata = new("AnnotatedDataFrame", data = scanDates, 
                varMetadata = data.frame(labelDescription = colnames(scanDates), 
                  row.names = colnames(scanDates)))
            protocolData(obj) = protocoldata
        }
        message("Determining chip type from IDAT protocolData...")
        if (ChipType == "BeadChip 12x1") {
            annotation(obj) = "IlluminaHumanMethylation27k"
        }
        else if (ChipType == "BeadChip 12x8") {
            annotation(obj) = "IlluminaHumanMethylation450k"
        }
        else if (ChipType == "BeadChip 8x5") {
            annotation(obj) = "IlluminaHumanMethylationEpic"
        }
    }
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
