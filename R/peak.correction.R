fuks <-
function (data, anno) {
   ds <- grep( 'DESIGN', colnames(anno) )
   stopifnot ( length(ds) == 1 )
   TI  <- anno[, ds] == "I"
   TII <- anno[, ds] == "II"
   corrected.data <- apply(data, 2, function(B) {
      C <- B ##LS
      SI <- summits(B[TI])
      SII <- summits(B[TII])
      C[TI] <- correctI(as.vector(B[TI]), SI, SII)
      C[TII] <- correctII(as.vector(B[TII]), SI, SII)
      return(C)
    })
    #row.names(corrected.data) <- c(row.names(data[TI, ]), row.names(data[TII, 
    #    ]))
    return(corrected.data)
}
