\name{BMIQ}
\alias{BMIQ}
\alias{CheckBMIQ}
\alias{BMIQ,MethyLumiSet-method}
\alias{BMIQ,MethylSet-method}
\alias{BMIQ-methods}
\alias{BMIQ,ANY-method}
\title{
Beta-Mixture Quantile (BMIQ) Normalisation method for Illumina 450k arrays
}
\description{
BMIQ is an intra-sample normalisation procedure, correcting the bias of type-2 probe values. BMIQ uses a 3-step procedure: (i) fitting of a 3-state beta mixture model, (ii) transformation of state-membership probabilities of type2 probes into quantiles of the type1 distribution, and (iii) a conformal transformation for the hemi-methylated probes. Exact details can be found in the reference below.
}
\usage{
BMIQ(beta.v, design.v, nL = 3, doH = TRUE, nfit = 50000, th1.v = c(0.2, 0.75), th2.v = NULL, niter = 5, tol = 0.001, plots = TRUE, sampleID = 1, pri=TRUE)
\S4method{BMIQ}{MethyLumiSet}(beta.v, nL=3, doH=TRUE, nfit=5000, th1.v=c(0.2,0.75), th2.v=NULL, niter=5, tol=0.001, plots=FALSE,  pri=FALSE ) 
CheckBMIQ(beta.v, design.v, pnbeta.v)
}
\arguments{
  \item{beta.v}{vector consisting of beta-values for a given sample, or a MethyLumiSet or MethylSet containing multiple samples. For the MethyLumiSet and MethylSet methods, this is the only required argument, and the function will be run on each sample.}
  \item{design.v}{corresponding vector specifying probe design type (1=type1,2=type2). This must be of the same length as beta.v and in the same order.}
  \item{nL}{number of states in beta mixture model. 3 by default. At present BMIQ only works for nL=3.}
  \item{doH}{perform normalisation for hemimethylated type2 probes. These are normalised using an empirical conformal transformation and also includes the left-tailed type2 methylated probes since these are not well described by a beta distribution. By default TRUE.}
  \item{nfit}{number of probes of a given design type to use for the fitting. Default is 50000. Smaller values (~10000) will make BMIQ run faster at the expense of a small loss in accuracy. For most applications, 5000 or 10000 is ok.}
  \item{th1.v}{thresholds used for the initialisation of the EM-algorithm, they should represent buest guesses for calling type1 probes hemi-methylated and methylated, and will be refined by the EM algorithm. Default values work well in most cases.}
  \item{th2.v}{thresholds used for the initialisation of the EM-algorithm, they should represent buest guesses for calling type2 probes hemi-methylated and methylated, and will be refined by the EM algorithm. By default this is null, and the thresholds are estimated based on th1.v and a modified PBC correction method.}
  \item{niter}{maximum number of EM iterations to do. This number should be large enough to yield good fits to the type1 distribution. By default 5.}
  \item{tol}{tolerance convergence threshold for EM algorithm. By default 0.001.}
  \item{plots}{logical specifying whether to plot the fits and normalised profiles out. By default TRUE.}
  \item{sampleID}{the ID of the sample being normalised.}
  \item{pri}{logical: print verbose progress information?}
  \item{pnbeta.v}{BMIQ normalised profile.}
}
\details{
Full details can be found in the reference below. Note: these functions require the RPMM package, not 
currently a dependency of the wateRmelon package.
}
\value{
Default method: A list with following entries:
\item{nbeta}{the normalised beta-profile for the sample}
\item{class1}{the assigned methylation state of type1 probes}
\item{class2}{the assigned methylation state of type2 probes}
\item{av1}{the mean beta-values for the nL states for type1 probes}
\item{av2}{the mean beta-values for the nL states for type2 probes}
\item{hf}{the estimated "Hubble" dilation factor used in the normalisation of hemi-methylated probes}
\item{th1}{estimated thresholds for calling unmethylated and methylated type1 probes}
\item{th2}{estimated thresholds for calling unmethylated and methylated type2 probes}

MethyLumiSet method: A methyLumiSet object

}
\references{
Teschendorff AE, Marabita F, Lechner M, Bartlett T, Tegner J, Gomez-Cabrero D,
Beck S. A Beta-Mixture Quantile Normalisation method for correcting probe design 
bias in Illumina Infinium 450k DNA methylation data. Bioinformatics. 2012 Nov 21.
}
\author{
Andrew Teschendorff, MethyLumiSet method by Leo Schalkwyk Leonard.Schalkwyk@kcl.ac.uk
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# library(RPMM)
# data(melon)
# BMIQ(melon,nfit=100)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
