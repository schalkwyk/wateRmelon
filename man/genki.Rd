\name{genki}
\alias{genki}
\title{
SNP derived performance metrics for Illumina 450K DNA methylation arrays.
}
\description{
A very simple genotype calling by one-dimensional K-means clustering is 
performed on each SNP, and for those SNPs where there are three genotypes, the squared deviations are summed for each genotype (similar 
to a standard deviation for each of allele A homozygote, heterozygote and 
allele B homozygote).  By default these are further divided by the square root
of the number of samples to get a standard error-like statistic.
}
\usage{
genki(bn, g = getsnp(rownames(bn)), se = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bn}{
a matrix of beta values(default method), a \code{MethyLumiSet} object (\code{methylumi} package), a \code{MethylSet} or \code{RGChannelSet} object (\code{minfi} package) or a \code{exprmethy450} object (\code{IMA} package).
}
  \item{g}{
vector of SNP names
}
  \item{se}{
TRUE or FALSE specifies whether to calculate the standard error-like statistic
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
There are 65 well-behaved SNP genotyping probes included on the array. 
These each produce a distribution of betas with tight peaks for the three 
possible genotypes, which will be broadened by technical variation between
samples.  The spread of the peaks is thus usable as a performance metric.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a vector of 3 values for the dispersion of the three genotype
peaks (AA, AB, BB : low, medium and high beta values)

}
\references{
Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
}
\author{
Leonard.Schalkwyk@kcl.ac.uk
}

\note{
Corrected \code{RGChannelSet} methods - 12/10/2015
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  #MethyLumiSet method
     data(melon)
     genki(melon)

  #MethyLumiSet method after normalization
     melon.dasen <- dasen(melon)
     genki(melon.dasen)
}

