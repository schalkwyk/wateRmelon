fuks <-
function (data, anno) {
   ds <- grep( 'DESIGN', colnames(anno) )
   stopifnot ( length(ds) == 1 )
   TI  <- anno[, ds] == "I"
   TII <- anno[, ds] == "II"
   corrected.data <- apply(data, 2, function(B) {
      SI <- summits(B[TI])
      SII <- summits(B[TII])
      BI <- correctI(as.vector(B[TI]), SI, SII)
      BII <- correctII(as.vector(B[TII]), SI, SII)
      return(c(BI, BII))
    })
    row.names(corrected.data) <- c(row.names(data[TI, ]), row.names(data[TII, 
        ]))
    return(corrected.data)
}
