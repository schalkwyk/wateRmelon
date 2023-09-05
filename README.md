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

You can read in EPICv2 array IDAT files, but I am curently testing and working through a number 
of functions that are broken by the changed probe naming convention.

## 1. Installation


**Install from Github**
```R
## Make sure 'devtools' is installed in your R
# install.packages("devtools")
devtools::install_github("schalkwyk/wateRmelon")
```
