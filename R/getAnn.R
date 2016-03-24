


# get the current Illumina annotation file 


aoget <- function(url= paste(
   'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/',
   'HumanMethylation450_15017482_v1-2.csv', sep='')
   ) {
   op <- options(stringsAsFactors=FALSE)
   ao <- read.csv(url, skip =7 )
   options(op)
   
}
