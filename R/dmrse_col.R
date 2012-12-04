dmrse_col <-
function(betas, idmr=iDMR()) {  # formerly SDR - between probes
    idmr <- idmr[idmr %in% rownames(betas)]
    sd(rowMeans(betas[idmr,], na.rm=T)) / sqrt(dim(betas)[2])
}
