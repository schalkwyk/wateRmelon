# Create Arbitrary Manifest using minfi engine:
# This needs to become much much more sophisticated...
.createAnnotation <- function(object){
	rn = rownames(object)
	message('Attempting to guess annotation, usually will only work on full-unfiltered dataset!')
	mode <- ifelse(length(rn) < 500000, "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19") # Weird force to handle RGChannelSets... Why would you use this on an RGChannelSet when you have get annotation. This function is only really intended for use on MethyllumiSets
	annoObj <-  minfi::getAnnotationObject(mode)
	all_data <- minfi:::.availableAnnotation(annoObj)$defaults
	new_data <- do.call(cbind, lapply(all_data, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
	}))
	new_data <- new_data[rn, ] # SNP probes will be missing, and be NA’d
	rownames(new_data) <- rn
	return(new_data)
}