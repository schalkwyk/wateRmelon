#New BScon_minfi after modifications of 17.09.2015
## BScon function adapted to minfi package
bscon_minfi <- function(RGsetEx){

#getting the values from the green channel
(csp.green <- function (RGsetEx, controls = c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II"))
{
  minfi:::.isRGOrStop(RGsetEx)
  r <- getRed(RGsetEx)
  g <- getGreen(RGsetEx)
  sapply (controls, function( controlType ) {

   ctrlAddress <- try (getControlAddress(RGsetEx, controlType = controlType), silent = T)
   if (!inherits (ctrlAddress, 'try-error')){ctrlAddress <- getControlAddress(RGsetEx, controlType = controlType)}
   else
     stop ("450k QC data could not be found")


   g[ctrlAddress, ]
  })})

green <- csp.green(RGsetEx)

#Getting values from the red channel
(csp.red <- function (RGsetEx, controls = c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II"))
{
  minfi:::.isRGOrStop(RGsetEx)
  r <- getRed(RGsetEx)
  g <- getGreen(RGsetEx)
  sapply (controls, function( controlType ) {

    ctrlAddress <- getControlAddress(RGsetEx, controlType = controlType)

    r[ctrlAddress, ]
  })})

red <- csp.red (RGsetEx)

#selecting only the Bisulfite conversion I values from both green and red
bsI.green <- green$`BISULFITE CONVERSION I`
bsI.red <- red$`BISULFITE CONVERSION I`
#selecting only the Bisulfite conversion II values from both green and red
bsII.green <- green$`BISULFITE CONVERSION II`
bsII.red <- red$`BISULFITE CONVERSION II`

# calculate BS conv type I betas as an example of using an index vector
if(nrow(bsI.green) > 11){ # 450K
  BSI.betas <- rbind(bsI.green[1:3,], bsI.red[7:9,])/((rbind(bsI.green[1:3,], bsI.red[7:9,])) + rbind(bsI.green[4:6,], bsI.red[10:12,]))
} else { # EPIC
  BSI.betas <- rbind(bsI.green[1:2,], bsI.red[6:7,])/((rbind(bsI.green[1:2,], bsI.red[6:7,])) + rbind(bsI.green[3:4,], bsI.red[ 8:9 ,]))
}

#calculation of BS con in Type II data
BSII.betas <- bsII.red/(bsII.red + bsII.green)

apply(rbind(BSI.betas, BSII.betas), 2, median)*100 ## this is the value you are interested in
}
