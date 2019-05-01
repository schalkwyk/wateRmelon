\name{betaqn-exprmethy450-methods}

\alias{betaqn,exprmethy450-method}
\alias{fuks,exprmethy450-method}

\title{
Calculate normalized betas from exprmethy450 of Illumina 450K methylation arrays
}  
\description{ 
Quantile normalize betas from exprmethy450 objects 
}

%\section{methods}{
%\describe{
%betaqn( bn ) 
%fuks( bn ) 
%}
%}


\arguments{

  \item{bn}{
An exprmethy450 object.
   }

  \item{fudge}{
value added to total intensity to prevent denominators close to zero when
calculating betas
   }
}
\details{
\bold{betaqn}	quantile normalizes betas 
}
\value{
exprmethy450 object of the same shape and order as bn. 

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
2012, 4:325-341

[4] Maksimovic J, Gordon L, Oshlack A: SWAN: Subset quantile 
Within-Array Normalization for Illumina Infinium HumanMethylation450 
BeadChips. Genome biology 2012, 13(6):R44

}
\author{
Leonard.Schalkwyk@kcl.ac.uk
}
%%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
