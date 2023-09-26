# Create Arbitrary Manifest using minfi engine: This needs to become much much
# more sophisticated...

.createAnnotation <- function(object) {
    rn = rownames(object)
    message("guessing array type from number of features ...")
    mode <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
    if(length(rn) < 5e+05) mode <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    if(length(rn) > 9e+05) mode <- "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
    message("using ", mode)
    annoObj <- minfi::getAnnotationObject(mode)
    all_data <- minfi:::.availableAnnotation(annoObj)$defaults
    new_data <- do.call(cbind, lapply(all_data, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
    }))
    new_data <- new_data[rn, ]  # SNP probes will be missing, and be NA’d
    rownames(new_data) <- rn
    return(new_data)
}



