#beadcount function, creates matrix with NAs representing probes with beadcount <3 from Extended RG Channel Set
beadcount<-function(x){
	#select out bead count dataframe
	assayDataElement(x, "NBeads")->nb

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
    	bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$Address,]
       
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


beadc<-function(x){
	length(which(is.na(x)=="TRUE"))}


# modified pfilter that only requires pn and bc, and applies filtering to any matrix given
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
