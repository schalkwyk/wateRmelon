# wateRmelon
 Illumina 450 and EPIC methylation array normalization and metrics

## Software status

| Resource:     | Bioconductor  (Release)      | Bioconductor (Devel)    |
| ------------- | ------------------- | ------------- |
| _Platforms:_  | _Multiple_          | _Multiple_    |
| R CMD check   | <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/wateRmelon/"><img border="0" src="http://bioconductor.org/shields/build/release/bioc/wateRmelon.svg" alt="Build status"></a></br>|<a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/wateRmelon/"><img border="0" src="http://bioconductor.org/shields/build/devel/bioc/wateRmelon.svg" alt="Build status"></a>

Currently on Github I'm doing some maintenance, in particular in order to work with EPICv2 arrays.
If you install:

   https://github.com/jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38

and

  https://github.com/jokergoo/IlluminaHumanMethylationEPICv2manifest

....  as well as this version of wateRmelon, you can read in EPICv2 array IDAT files with wateRmelon and bigmelon.

Most functions including dasen and the other Pidsley et al 2013 normalisers, genki, dmrse,
 and agep work at least in wateRmelon but have only been superficially tested, additional testingand improvement is underway.  Bug reports are welcome at lschal@essex.ac.uk or via github.


## 1. Installation


**Install from Github**
```R
## Make sure 'devtools' is installed in your R
# install.packages("devtools")
devtools::install_github("schalkwyk/wateRmelon")
```
