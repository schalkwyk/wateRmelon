# Leo 26 Nov 2012
# wrapper for combine method

combo <- function (...){
   if(!library(methylumi, logical.return=TRUE, quietly=TRUE)){
      stop('can\'t load methylumi package')
   }
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

