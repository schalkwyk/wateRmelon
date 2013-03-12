\name{seabird}
\alias{seabird}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate ROC area-under-curve for X-chromosome sex differences (internal function for calculating the seabi metric)
}
\description{
This is a wrapper for the prediction and performance functions from
the ROCR package that takes a vector of p-values and a vector of true
or false for being on the X.  See \code{seabi} function which does everything.
}
\usage{
seabird(pr, stop = 1, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pr}{
a vector of p-values, such as calculated by \code{seabird}
}
  \item{stop}{
fraction for partial area under curve.  For example 0.1 gives you the
area for the lowest 10\% of p-values. 
}
  \item{X}{
logical vector the same length as pv, true for features mapped to X-chromosome
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns an area value between 0 and 1, where 1 is the best possible performance. 
}
\references{
Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: 
A data-driven approach to preprocessing Illumina 450K methylation 
array data (submitted)
%% ~put references to the literature/web site here ~
}
\author{
Leonard C Schalkwyk 2012
Leonard.Schalkwyk@kcl.ac.uk
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
