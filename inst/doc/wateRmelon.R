### R code from vignette source 'wateRmelon.Rnw'

###################################################
### code chunk number 1: UnevaluatedCode (eval = FALSE)
###################################################
## install.packages('ROCR', 'matrixStats')
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install( 'limma', 'minfi',
##    'IlluminaHumanMethylation450kmanifest',
##    'methylumi', 'lumi')


###################################################
### code chunk number 2: UnevaluatedCode (eval = FALSE)
###################################################
## install.packages('wateRmelon_0.9.9.tar.gz', repos=NULL, type='source')


###################################################
### code chunk number 3: code-block
###################################################
library('wateRmelon')
# load in melon dataset
data (melon)
# display dimensions of data matrix
dim(melon)
# quality filter using default thresholds
melon.pf<-pfilter(melon)
# preprocess using our best method
melon.dasen.pf <- dasen(melon.pf)


###################################################
### code chunk number 4: dmrse
###################################################
# calculate iDMR metrics on QC'd betas
dmrse_row(melon.pf)
# calculate iDMR metrics on QC'd and preprocessed betas
dmrse_row(melon.dasen.pf)
# slightly lower (better) standard errors


###################################################
### code chunk number 5: genki
###################################################
# calculate SNP metrics on QC'd betas
genki(melon.pf)
# calculate SNP metrics on QC'd and preprocessed betas
genki(melon.dasen.pf)
# slightly lower (better) standard errors


###################################################
### code chunk number 6: seabi
###################################################
# calculate X-chromosome metrics on QC'd betas
seabi(melon.pf, sex=pData(melon.pf)$sex, X=fData(melon.pf)$CHR=='X')
# calculate X-chromosome metrics on QC'd and preprocessed betas
seabi(melon.dasen.pf,
   sex=pData(melon.dasen.pf)$sex,
   X=fData(melon.dasen.pf)$CHR=='X'
)


###################################################
### code chunk number 7: UnevaluatedCode (eval = FALSE)
###################################################
## library(methylumi)
## melon <- methyLumiR('finalreport.txt')


###################################################
### code chunk number 8: UnevaluatedCodeecc (eval = FALSE)
###################################################
## mlumi <- readEPIC('path/to/directory')


###################################################
### code chunk number 9: IncludeGraphic
###################################################
boxplot(log(methylated(melon)), las=2, cex.axis=0.8 )


###################################################
### code chunk number 10: IncludeGraphic
###################################################
boxplot(log(unmethylated(melon)), las=2, cex.axis=0.8 )


###################################################
### code chunk number 11: outlyx
###################################################
outlyx(melon) # Can take some time on large data-sets


###################################################
### code chunk number 12: bscon
###################################################
bsc <- bscon(melon)
hist(bsc)


###################################################
### code chunk number 13: pfilter
###################################################
melon.pf <- pfilter(melon)


###################################################
### code chunk number 14: dasen
###################################################
melon.dasen.pf <- dasen(melon.pf)


###################################################
### code chunk number 15: agep
###################################################
data(coef)
agep(melon.dasen.pf, coeff= coef, method='horvath')


###################################################
### code chunk number 16: UnevaluatedCodeecc (eval = FALSE)
###################################################
## # Code will not work with melon as melon only has a subset of probes
## estimateCellCounts.wmln(melon.dasen.pf)


###################################################
### code chunk number 17: qual
###################################################
normv <- qual(betas(melon.dasen.pf), betas(melon.pf))
plot(normv[,1:2])


###################################################
### code chunk number 18: pwod
###################################################
melon.pwod.dasen.pf <- pwod(betas(melon.dasen.pf))


###################################################
### code chunk number 19: workflow
###################################################
data(melon)
# load in melon dataset
# Optional read-in data using:
# melon <- readEPIC('path/to/idats')

outliers <- outlyx(melon, plot=F)
sum(outliers$out)
sum(bscon(melon)<85)
# Check for outliers

melon.pf<-pfilter(melon)
# perform QC on raw data matrix using default thresholds
melon.dasen.pf<-dasen(melon.pf)
# preprocess using our best method

qual <- qual(betas(melon.dasen.pf), betas(melon.pf))
# Check for bad samples (again)

sex  <- pData(melon.dasen.pf)$sex
# extract phenotypic information for test
bet<-betas(melon.dasen.pf)
# extract processed beta values
melon.sextest<-sextest(bet,sex)
# run t-test to idenitify sex difference

agep(melon.dasen.pf)
# Check ages (see if they match up e.t.c)

melon.pwod <- pwod(bet)
# Clean up data for statistical testing


