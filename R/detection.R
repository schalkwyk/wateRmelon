#' detectP -  minfi-style detection p-values from methylumi objects
# the original detection pval was produced by Illumina's software and 
# afaik not documented.  It's basically a way of thresholding signal
# intensity to filter and/or check for numbers of failed probes and
# its use is enshrined in pfilter().

# There are several confusingly named functions in methylumi and minfi
# including accessors of detection pvalues.  

# minfi's detectionP calculates a reasonable t-test like detection
# pvalue using the mean (actually median) and sd (actually mad) of 
# the negative control probes as a background value.  This only has 
# methods for RGset (ie unpreprocessed) objects.  They calculate these 
# stats for red and green and then double these to test against typeI
# and add red and green for type 2.  Here I'll just use m+u both sides.
# This may behave slightly different, I'll test afterwards.
#
# note  that there are some functions including readEPIC that just 
# write betas into the pvals slot as an expedient placeholder.



detectP <- function(mlo){

   qc <- QCdata(mlo) 
   neg <- featureData(qc)$Type=='NEGATIVE'
   
   # mean and mad controls for each sample 
   
   control <- methylated(qc)[neg,,drop=FALSE] + 
            unmethylated(qc)[neg,,drop=FALSE]
   Mean <- colMedians(
      x=control, rows=NULL, cols=NULL, na.rm=T, use.names=TRUE 
   )
   Mad  <- colMads(
      x=control, rows=NULL, cols=NULL, na.rm=T, use.names=TRUE 
   )

   # test.  Note 'good' samples have low numbers, eg p < .05.
   # also note the transpose is needed because of recycling:
   # we're giving pnorm a row of means and mads.  Minfi instead 
   # loops through the columns.
   
   ints <- methylated(mlo) + unmethylated(mlo)
   t(pnorm( q=t(ints), mean=Mean, sd=Mad, lower.tail=FALSE, log.p=FALSE))

}
