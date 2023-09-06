# Leo 26 Nov 2012
# wrapper for combine method



#' Combine MethyLumiSet objects
#' 
#' This is a wrapper for combining different MethyLumiSet objects.
#' 
#' This is a wrapper for \code{methylumi::combine}, which works around a name
#' clash with a different combine function from the \code{gdata} package, and
#' also a bug in \code{methylumi::combine}.
#' 
#' @param \dots Eventually, any number of MethyLumiSet objects. Currently only
#' guaranteed for 2 objects.
#' @return a \code{MethyLumiSet}.  The \code{assayData}, \code{QCdata},
#' \code{experimentData}, \code{protocolData} and \code{phenoData} are joined
#' on \code{sampleName} .  \code{featureData} and annotation are taken from the
#' object given in the first argument
#' @note the function uses \code{sampleNames} and gets rid of duplicates.
#' Numeric sampleNames cause problems (and are a Bad Idea anyway).  They should
#' be turned into names with \code{make.names()} first.
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>
#' @seealso \code{\link{as.methylumi}}
#' @references [1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk
#' LC: A data-driven approach to preprocessing Illumina 450K methylation array
#' data (submitted)
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' library(wateRmelon)
#' data(melon)
#' ## pretend we have two different data sets
#' melon
#' pelon <- melon
#' sampleNames(pelon) <- gsub('^6', 7, sampleNames(pelon))
#' combo(melon, pelon)
#' 
#' 
#' 
#' @export combo
combo <- function (...){
   history.submitted <- as.character(Sys.time())   
   vict  <- list(...)
   asda  <- lapply( vict, assayData )
   peda  <- lapply( vict, pData )
   prda  <- lapply( vict, protocolData )
   exda  <- lapply( vict, experimentData )
   qcda  <- lapply( vict, function(x) assayData(QCdata(x)))
   clod  <- new('MethyLumiSet')
   assayData(clod)      <- do.call(methylumi::combine,asda) 
   pData(clod)          <- do.call(methylumi::combine,peda) 
   protocolData(clod)   <- do.call(methylumi::combine,prda) 
   experimentData(clod) <- do.call(methylumi::combine,exda) 
    QCdata(clod) <- new("MethyLumiQC", assayData = do.call(methylumi::combine, qcda))
   fData(clod)     <- fData(vict[[1]] )
   annotation(clod)<- annotation(vict[[1]] )
   history.finished <- as.character(Sys.time())
   clod@history<- data.frame(
      submitted = history.submitted, 
      finished  = history.finished, 
      command   = paste(as.character(sys.call()), collapse=' ')
   )
   clod
}

