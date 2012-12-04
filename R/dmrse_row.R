dmrse_row <-
function(betas, idmr=iDMR()) {  # formerly SDC - between arrays
    idmr <- idmr[idmr %in% rownames(betas)]
    message ( paste ( length(idmr) , 'iDMR data rows found\n' ))
    sd(colMeans(betas[idmr,], na.rm=T)) / sqrt(dim(betas)[2])
}
