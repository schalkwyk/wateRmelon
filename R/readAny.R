#' readPepo - read (any kind of) Illumina DNA methylation array idat 
#' files into a methylumi object 
#'
#' Pepo is a botanical term for any melon-like fruit. Given the 
#' appropriate manifest file, this function should be able to read
#' any of the Illumina Infinium DNAm arrays including 450k, EPIC, EPIC2,
#' MSA and Mouse.
#'
#' @param idatdir the directory with the idatfiles. Currently only handle one directory. 
#' @param filelist optional list of idat files to process.
#' @param barcodelist optional list of barcodes to process.
#' @param manifest name of a IlluminaMethylationManifest object or a csv format manifest.  If missing,  will run idet() on one of the idat files.
#' @param parallel try to use multiple cores.
#' @param n keep beadcounts.
#' @param keep out-of-band (OOB) or opposite-channel signals
#' @param pdat optional data.frame describing the samples.
#' @param two are there two different assay types (true of human methylation arrays except 27k)
#' @return  A ‘MethyLumiSet’ object.
#' @examples
#'
#' @export


readPepo <- function (  idatdir='.' , filelist=NULL,
                        barcodelist=NULL, manifest=NULL, 
                        parallel=F, n=F, pdat=NULL, oob=F,
                        two=TRUE
                     ){


# check for matched red and green files

   if (!is.null(barcodelist)){
      if (!is.null(filelist)) {
         message("Both file and barcode lists specified, using barcodes")
      }   
      filelist <- c(
         paste(idatdir, '/', barcodelist, '_Red.idat', sep=''), 
         paste(idatdir, '/', barcodelist, '_Grn.idat', sep='')
      )
   } else {
      if (is.null(filelist)) filelist <- dir(idatdir, patt='idat') 
      # this should use basename
      barcodelist <- unique(substr(filelist,0,nchar(filelist)-9))  
      # this should prepend idatdir
      filelist <- c(
         paste(idatdir, '/', barcodelist, '_Red.idat', sep=''), 
         paste(idatdir, '/', barcodelist, '_Grn.idat', sep='')
      )
   }

# nothing:  list of idat files present and barcodes extracted from them
# barcodes or both: file list expanded from barcode list
# file list: given file list and barcodes from it, re-expand to catch unpaired 

   check <- sapply(filelist, file.exists)
   if (any(!check)){ stop ( "file(s)", filelist[!check], "missing" )}

# identify and if necessary process manifest.  
# manifest should be a IlluminaMethylationManifest obj when done.

   # name object from loaded package
   if(is.null(manifest)){ manifest <- idet(filelist[1])[3] } 
   # file or name of object from loaded package
   if(is.character(manifest)){   
      if(!is.na(file.info(manifest)$size)){ manifest <- canno(manifest) } 
      else{ manifest <- eval(parse(text=manifest)) }
   }
   stopifnot( is(manifest, 'IlluminaMethylationManifest'))  
   # manifest now variable in function env; stash copy in global
   .manifest <<- manifest

# now we do the thing (copied from methylumIDATepic)

   mats  <- IDATsToMatrices2(barcodelist, parallel = parallel, idatPath = idatdir)
   dats  <- DataToNChannelSet2(mats, IDAT = T, parallel = parallel, force = FALSE)
   annotation(dats)<- '.manifest'  # looks like correct NChannelSet
   mlumi <- NChannelSetToMethyLumiSet2(   
      dats, parallel = parallel, oob = oob, n = n
   )

   if(is.null(pdat)) {
      pdat           <- data.frame(barcode = as.character(barcodelist))
      rownames(pdat) <- pdat$barcode
      pData(mlumi)   <- pdat
   }  else {
        pData(mlumi) <- pdat
   }
   if (!is.null(mlumi@QC)) {
      sampleNames(mlumi@QC) = sampleNames(mlumi)}
   colnames(mlumi) <- as.character(colnames(mlumi))
   
#    return(mlumi[sort(featureNames(mlumi)), ])
#  we need to tag the NChannelSet with a meaningful Annotation tag that points to the 
#  appropriate manifest.  minfi -> readEPIC does this indirectly with a lookup table 
#  implmented as 3 switch statements in getMethylationBeadMappers2() and another in  
#  generateManifest(), which is called by getMethylationBeadMappers2().  We want to 
#  specify the manifest object directly.  
#  Ideally we would also clean up  getMethylationBeadMappers2() 
#  for the moment I am going to use the expedient of setting the Annotation to ".manifest"
#  and <<- the manifest to a hidden var .manifest in globalenv

   annotation(mlumi) <- '.manifest'
   mlumi




#echo top block of manifest
#
#IDATsToMatrices2 
#DataToNChannelSet2 (IDAT=F) #should make NChannelSet without setting annotation
#==
#option A:  
#process manifest file into object required by mergeProbeDesigns2 (and optionally store as rda ) and use NChannelSetToMethyLumiSet2
#
#option  B:
#re-implement processing
#==

}
