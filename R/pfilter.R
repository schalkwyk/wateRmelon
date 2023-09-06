#beadcount function, creates matrix with NAs representing probes with beadcount <3 from Extended RG Channel Set


#' Creates matrix of beacounts from minfi data.
#' 
#' Creates matrix of beacounts from data read in using the minfi package.  NAs
#' represent probes with beadcount <3.  An Extended RG Channel Set is required
#' for this function to work.
#' 
#' 
#' @param x 450K methylation data read in using minfi to create an Extended RG
#' Channel Set
#' @return A matrix of bead counts with bead counts <3 represented by NA for
#' use in the pfilter function for quality control
#' @note The beadcount function is internal to the pfilter function
#' @author Ruth.Pidsley
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>
#' @references [1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk
#' LC: A data-driven approach to preprocessing Illumina 450K methylation array
#' data (submitted)
#' @export beadcount
beadcount<-function(x){
	#select out bead count dataframe
    getNBeads(x) -> nb
	#match rownames of beadcount dataframe to addresses
	getProbeInfo(x,type="I")->typeIadd
	match(typeIadd$AddressA,rownames(nb))->typeImatchA
	match(typeIadd$AddressB,rownames(nb))->typeImatchB

	#match rownames of beadcount dataframe to addresses
	getProbeInfo(x,type="II")->typeIIadd
	match(typeIIadd$Address,rownames(nb))->typeIImatch

	nb->nbcg

    	locusNames <- getManifestInfo(x, "locusNames")
    	bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
        dimnames = list(locusNames, sampleNames(x)))
    
    	TypeII.Name <- getProbeInfo(x, type = "II")$Name
    	bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$AddressA,]
       
	TypeI <- getProbeInfo(x, type = "I")

    	bc_temp->bcB
    	bc_temp->bcA		
    
    	bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB,]
    	bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA,]

	which(bcB<3)->bcB3
	which(bcA<3)->bcA3
    	bcA->bcA2
    	bcB->bcB2

    	bcA2[bcA3]<-NA
    	bcB2[bcB3]<-NA
    	bcA2[bcB3]<-NA
    	bcB2[bcA3]<-NA
	
    	data.frame(bcA2)->bcM
    	data.frame(bcA2)->bcU
    	data.frame(bcA2)->bc
    	bc
}




#' Calculates the number of samples with bead count <3 for each probe in matrix
#' of bead count values
#' 
#' Calculates the number of samples with bead count <3 for each probe in matrix
#' of bead count values.
#' 
#' 
#' @param x matrix of bead count values returned by the beadcount function
#' @return Vector of number of samples with bead count <3 for each probe
#' @note The beadc function is internal to the pfilter function
#' @author Ruth.Pidsley
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>
#' @references [1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk
#' LC: A data-driven approach to preprocessing Illumina 450K methylation array
#' data (submitted)
#' @export beadc
beadc<-function(x){
	length(which(is.na(x)=="TRUE"))}


# modified pfilter that only requires pn and bc, and applies filtering to any matrix given


#' Basic data filtering for Illumina 450 methylation data
#' 
#' The pfilter function filters data sets based on bead count and detection
#' p-values.  The user can set their own thresholds or use the default pfilter
#' settings.  pfilter will take data matrices of beta values, signal
#' intensities and annotation data, but will also take methylumi (MethyLumiSet)
#' or minfi (RGChannelSetExtended) objects.  However it has come to our
#' attention that data read in using the various packages and input methods
#' will give subtly variable data output as they calculate detection p-value
#' and beta values differently, and do/don?t give information about beadcount.
#' The pfilter function does not correct for this, but simply uses the
#' detection p-value and bead count provided by each package.
#' 
#' 
#' @aliases pfilter-methods pfilter,MethyLumiSet-method
#' pfilter,RGChannelSetExtended-method pfilter
#' @param mn matrix of methylated signal intensities, each column representing
#' a sample (default method), or an object for which a method is available e.g
#' MethyLumiSet or RGChannelSetExtended.  N.B. Bead count filtering will not
#' work unless data read in as an RGGhannelSetExtended rather than an
#' RGChannelSet.
#' @param un matrix of unmethylated signal intensities, each column
#' representing a sample (default method) or NULL when mn is a MethyLumiSet or
#' RGChannelSetExtended object
#' @param bn matrix of precalculated betas, each column representing a sample,
#' or NULL when mn is a MethyLumiSet or RGChannelSetExtended object
#' @param da annotation data frame, such as x@featureData@data #methylumi
#' package, or NULL when mn is a MethyLumiSet or RGChannelSetExtended object
#' @param pn matrix of detection p-values, each column representing a sample, a
#' MethyLumiSet or RGChannelSetExtended object
#' @param bc matrix of arbitrary values, each column representing a sample and
#' eeach row representing a probe, in which "NA" represents beadcount <3, or
#' NULL when mn is a MethyLumiSet or RGChannelSetExtended object
#' @param perCount remove sites having this percentage of samples with a
#' beadcount <3, default set to 5
#' @param pnthresh cutoff for detection p-value, default set to 0.05
#' @param perc remove samples having this percentage of sites with a detection
#' p-value greater than pnthresh, default set to 1
#' @param pthresh remove sites having this percentage of samples with a
#' detection p-value greater than pnthresh, default set to 1
#' @param logical.return If it is TRUE, FALSE or TRUE is returned to indicate
#' success
#' @return a filtered MethyLumiSet or RGChannelSetExtended object or
#' 
#' a list of the filtered matrices:
#' 
#' mn : methylated intensities
#' 
#' un : unmethylated intensities
#' 
#' bn : betas
#' 
#' da : feature data
#' @section Methods: \describe{ \item{list("signature(mn =
#' \"MethyLumiSet\")")}{ This is used for performing the pfilter method on
#' MethyLumiSet objects produced by methylumiR.  } \item{list("signature(mn =
#' \"RGChannelSetExtended\")")}{ This is used for performing the pfilter method
#' on RGChannelSetExtended objects produced by minfi.  } }
#' @author Ruth.Pidsley
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>
#' @references [1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk
#' LC: A data-driven approach to preprocessing Illumina 450K methylation array
#' data (submitted)
#' @examples
#' 
#' # MethyLumiSet method
#' data(melon)
#' melon.pf <- pfilter(melon)
#' 
#' @export pfilter
pfilter<-function (mn=NULL, un=NULL, bn=NULL, da=NULL, pn, bc, perCount = NULL, pnthresh = NULL, perc = NULL, pthresh = NULL, logical.return=FALSE){
		
	if (!is.null(list(pn, bc))) {

		if (is.null(perCount)) {
             		perCount = 5
             	}
      		
		if (is.null(pnthresh)) {
             		pnthresh = 0.05
             	}
             	
		if (!is.null(pnthresh)) {
                 	pnthresh = pnthresh
             	}
             	

		if (is.null(perc)) {
                 	perc = 1
          	}
             	
		if (!is.null(perc)) {
                 	perc = perc
             	}
             	
		if (is.null(pthresh)) {
                 	pthresh = 1
             	}
             	
		if (!is.null(pthresh)) {
                 	pthresh = pthresh
	        }
             	
		goodsamp <- (colSums(pn > pnthresh)) < ((nrow(pn) * perc)/100)
      
      bab <- apply(bc[,goodsamp],1,beadc)
      badbead <- which(bab>(ncol(bc) * perCount)/100)
      badbead_log<-bab>(ncol(bc) * perCount)/100

      bap  <- rowSums(pn[, goodsamp] > pnthresh) > ((ncol(pn) * pthresh/100))
      badp <- which(bap)
        	
		cat( length(which(goodsamp==FALSE)), "samples having", perc, "% of sites with a detection p-value greater than", pnthresh, "were removed","\n")
		cat("Samples removed:", names(which(goodsamp==FALSE)), "\n")
		cat(length(badbead), "sites were removed as beadcount <3 in", perCount, "% of samples","\n")
		cat(length(badp), "sites having", pthresh, "% of samples with a detection p-value greater than", pnthresh, "were removed", "\n")
		
      if(logical.return) return(list(probes=(!bap & !badbead_log), samples=goodsamp)) 

		if (!is.null(mn)) {
			mn <- mn[-c(badp, badbead), goodsamp]
		}
             	if (!is.null(un)) {
			un <- un[-c(badp, badbead), goodsamp]
		}
		if (!is.null(bn)) {
             		bn <- bn[-c(badp, badbead), goodsamp]
		}
		if (!is.null(da)) {
             		da <- da[-c(badp, badbead), ]
		}
		

     	        
	}
	list(mn = mn, un = un, bn = bn, da = da)
	
}
