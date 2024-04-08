###  LS Mar 2024 functions related to manifests and chip types  



#' read.manifest - read in csv format Illumina chip manifest files
#'
#' @param file 
#' @return  a list of of dataframes of data prepared for making IlluminaMethylationManifest
#' @details
#'   This function is probably not much use for calling directly, it mostly exists to be called by canno.
#'   

read.manifest <- function(file) {
    
   seps <- grep('^\\[', readLines(file) )
   manifest <- read.csv(file, skip=seps[2], nrows=seps[3]-seps[2]-2, as.is=T)
   manifest$Name ->  manifest$oldName  # with EPICv2 names are no longer unique!
   manifest$Name <- manifest$IlmnID
   TypeI <- manifest[
        manifest$Infinium_Design_Type == "I",
        c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel", "AlleleB_ProbeSeq" )]
   names(TypeI)[c(2, 3, 4)] <- c( "AddressA", "AddressB", "Color" )
   TypeI <- as(TypeI, "DataFrame")
   TypeI$ProbeSeqB <- DNAStringSet(TypeI$AlleleB_ProbeSeq)
   TypeI$nCpG <- as.integer(
        oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
   TypeI$nCpG[TypeI$nCpG < 0] <- 0L
       # there are 295 probes including 30 cpg probes and nv- and rs- annotated SNPs with no CG 

   TypeI$Name     <- as.character(TypeI$Name    ) 
   TypeI$AddressA <- as.character(TypeI$AddressA) 
   TypeI$AddressB <- as.character(TypeI$AddressB) 
   TypeI$Color    <- as.character(TypeI$Color   ) 
   TypeI$nCpG     <- as.integer  (TypeI$nCpG    )
    
   TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]   # unlike minfi we want these in the main data set

   TypeII <- manifest[
       manifest$Infinium_Design_Type == "II",
       c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
   names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
   TypeII <- as(TypeII, "DataFrame")
   TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
   TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
       # probeseqA don't contain CG because designed to anneal to BS DNA, complement of CG is RC
       # this correct if there are no other Rs in the sequences, haven't studied in detail, but 
       # there do seem to be more probes with Rs than RC. 
   TypeII$nCpG[TypeII$nCpG < 0] <- 0L

   TypeII$Name       <- as.character(TypeII$Name    ) 
   TypeII$AddressA   <- as.character(TypeII$AddressA)
   TypeII$nCpG       <- as.integer  (TypeII$nCpG    ) 
  
   TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]

    controls <- read.table(
        file = file,
        skip = seps[3],
        sep = ",",
        comment.char = "",
        quote = "",
        colClasses = c(rep("character", 5)))[, 1:5]
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    TypeControl <- as(TypeControl, "DataFrame")

   TypeControl$Type      <-as.character(TypeControl$Type   ) 
   TypeControl$Address   <-as.character(TypeControl$Address)


    list(
        manifestList = list(
            TypeI = TypeI,
            TypeII = TypeII,
            TypeControl = TypeControl,
            TypeSnpI = TypeSnpI,
            TypeSnpII = TypeSnpII),
        manifest = manifest,
        controls = controls)
}

#' canno  - process csv manifest into annotation object for illumina methylation preprocessing
#'
#' @param man name 
#' @return  IlluminaMethylationManifest
#' @details    
#' This is based on the scripts shipped with minfi annotation packages. It is based on existing 
#' Illumina Human Methylation csv format manifests, but because reverse-engineered, may require
#' updates to work on future products.
#'   
#' 

#


canno <- function (man= "EPIC-8v2-0_A1.csv", name=NULL){  

   maniTmp <-      read.manifest(man)
   manifestList <- maniTmp$manifestList
   if (is.null(name)) name <- basename(man)

  IlluminaMethylationManifest(
     TypeI = manifestList$TypeI,
     TypeII = manifestList$TypeII,
     TypeControl = manifestList$TypeControl,
     TypeSnpI = manifestList$TypeSnpI,
     TypeSnpII = manifestList$TypeSnpII,
     annotation = name)
  
}


#' idet - identify idats by a hash of the addresses
#'
#' @param idat cn be an idat file or the list produced by reading one with readIDAT()
#' @return  three strings: ChipType, an md5 hash of the MidBlock (address vector) and if known the annotation name
#' @details
#'   this function is a response to the fact that IlluminaHumanMethylationEPIC and EPICv2
#'   idats both have the ChipType "BeadChip 8x5" but different manifests.  They did have
#'   different numbers of addresses.  Subsequently we have had some confusion....
#'   This hash (certainly taken together with the ChipType) should be bombproof
#'   as an identifier.  Warning: this is slow.



idet <- function(idat){

   knowns <- c(
      `2f96172b55a56a146c47b4d48cdfb3e0`  = "IlluminaHumanMethylationEpic",    
      `d4e153ca79918396f3371bad7ef8082d`  = "IlluminaHumanMethylationEpic",    
      `e28439b7c9a50e3a0f11252f02d09c7e`  = "IlluminaHumanMethylationEpic",    
      `485216d14ac9190d13df19681ec775c7`  = "IlluminaHumanMethylationEpicv2"  

   )

   if(is.character(idat)) idat <- readIDAT(idat)
   write(idat$MidBlock, x<-R.utils::tmpfile())
   hesh <- tools::md5sum(x) ; unlink(x)
   names(hesh) <- NULL
   c(idat$ChipType, hesh, knowns[hesh])

}
