\name{smokp}
\alias{smokp}
\alias{smokp,MethylSet-method}
\alias{smokp,RGChannelSet-method}
\alias{smokp,MethyLumiSet-method}
\alias{smokp_cpgs}
\title{
  Smoking Prediction from methylomic expression data
}
\description{
  Predict smoking from samples using various methods
}
\usage{
  smokp(betas, method, sst)
}
\arguments{
  \item{betas}{
    Matrix of betas or MethyLumiSet or MethylSet object. Rows are Illumina IDs referring to CpG sites and Columns refer to samples or participants.
  }
  \item{method}{
    Currently: 'AHRR', 'McCartney', 'Maas', 'Sugden', 'Teschendorff', 'Yu', 'Gao', 'Yang', 'Zhang', 'Wen', 'Langdon', 'SSt', 'Packyears', 'Cessation', and 'All'. If "All" smokp will seek to predict smoking using all methods else will use the method specified. Default is "SSt". If 'Teschendorff', 'Yu', 'Gao', 'Yang' or 'Langdon' specified then smoking status is required.
  }
  \item{sst}{
    Named vector describing smoking status, coded as 'Current', 'Former', or 'Never', of participant for each sample (where names match rownames of betas).
  }}
  \value{
    Returns data frame of predicted smoking per sample.
  }
  \author{
    Original Functions: See References.
    
    wateRmelon: Tyler Gorrie-Stone, Leo Schalkwyk, Louis El Khoury
    
    smokp: Alexandria Andrayas
  }
  \references{
  Philibert, R.A., Beach, S.R. and Brody, G.H. Demethylation of the aryl hydrocarbon receptor repressor as a biomarker for nascent smokers. Epigenetics, 7:11, 2012, 1331-1338.
  
  Philibert, R., Hollenbeck, N., Andersen, E., McElroy, S., Wilson, S., Vercande, K., Beach, S.R., Osborn, T., Gerrard, M., Gibbons, F.X. and Wang, K. Reversion of AHRR demethylation is a quantitative biomarker of smoking cessation. Frontiers in psychiatry, 2016, 7, 55.
  
  Zeilinger, S., Kühnel, B., Klopp, N., Baurecht, H., Kleinschmidt, A., Gieger, C., Weidinger, S., Lattka, E., Adamski, J., Peters, A. and Strauch, K., Tobacco smoking leads to extensive genome-wide changes in DNA methylation. 2013, PloS one, 2013, 8:5, e63812.
  
  Elliott, H.R., Tillin, T., McArdle, W.L., Ho, K., Duggirala, A., Frayling, T.M., Smith, G.D., Hughes, A.D., Chaturvedi, N. and Relton, C.L. Differences in smoking associated DNA methylation patterns in South Asians and Europeans. Clinical epigenetics, 2014, 6:1, 1-10.
  
  Teschendorff, A.E., Yang, Z., Wong, A., Pipinikas, C.P., Jiao, Y., Jones, A., Anjum, S., Hardy, R., Salvesen, H.B., Thirlwell, C. and Janes, S.M. Correlation of smoking-associated DNA methylation changes in buccal cells with DNA methylation changes in epithelial cancer. JAMA oncology, 2015, 1:4, 476-485.
  
  Zhang, Y., Florath, I., Saum, K. U., & Brenner, H. Self-reported smoking, serum cotinine, and blood DNA methylation. Environmental research, 2016, 146, 395-403.
  
  Gao, X., Zhang, Y., Breitling, L.P. and Brenner, H. Relationship of tobacco smoking and smoking-related DNA methylation with epigenetic age acceleration. Oncotarget, 2016, 7:30, 46878.
  
  Zhang, Y., Schöttker, B., Florath, I., Stock, C., Butterbach, K., Holleczek, B. & Brenner, H. Smoking-associated DNA methylation biomarkers and their predictive value for all-cause and cardiovascular mortality. Environmental health perspectives, 2016, 124:1, 67-74.
  
  McCartney, D.L., Hillary, R.F., Stevenson, A.J., Ritchie, S.J., Walker, R.M., Zhang, Q., Morris, S.W., Bermingham, M.L., Campbell, A., Murray, A.D. and Whalley, H.C. Epigenetic prediction of complex traits and death. Genome biology, 2018, 19:1, 1-11. 
  
  Joehanes, R., Just, A.C., Marioni, R.E., Pilling, L.C., Reynolds, L.M., Mandaviya, P.R., Guan, W., Xu, T., Elks, C.E., Aslibekyan, S. and Moreno-Macias, H. Epigenetic signatures of cigarette smoking. Circulation: cardiovascular genetics, 2016, 9:5, 436-447.
  
  Sugden, K., Hannon, E.J., Arseneault, L., Belsky, D.W., Broadbent, J.M., Corcoran, D.L., Hancox, R.J., Houts, R.M., Moffitt, T.E., Poulton, R. and Prinz, J.A. Establishing a generalized polyepigenetic biomarker for tobacco smoking. Translational psychiatry, 2019, 9:1, 1-12.
  
  Gao, X., Jia, M., Zhang, Y., Breitling, L.P. and Brenner, H. DNA methylation changes of whole blood cells in response to active smoking exposure in adults: a systematic review of DNA methylation studies. Clinical epigenetics, 2015, 7, 1-10.
  
  Yu, H., Raut, J.R., Schöttker, B., Holleczek, B., Zhang, Y. and Brenner, H. Individual and joint contributions of genetic and methylation risk scores for enhancing lung cancer risk stratification: data from a population-based cohort in Germany. Clinical epigenetics, 2020, 12:1, 1-11.
  
  Yang, Y., Gao, X., Just, A.C., Colicino, E., Wang, C., Coull, B.A., Hou, L., Zheng, Y., Vokonas, P., Schwartz, J. and Baccarelli, A.A. Smoking-related DNA methylation is associated with DNA methylation phenotypic age acceleration: The veterans affairs normative aging study. International journal of environmental research and public health, 2019, 16:13, 2356.
  
  Bollepalli, S., Korhonen, T., Kaprio, J., Anders, S., & Ollikainen, M. EpiSmokEr: A robust classifier to determine smoking status from DNA methylation data. Epigenomics, 2019, 11:13, 1469-1486.
  
  Maas, S.C., Vidaki, A., Wilson, R., Teumer, A., Liu, F., van Meurs, J.B., Uitterlinden, A.G., Boomsma, D.I., de Geus, E.J., Willemsen, G. and van Dongen, J. Validated inference of smoking habits from blood with a finite DNA methylation marker set. European journal of epidemiology, 2019, 34, 1055-1074.
  
  Langdon, R.J., Yousefi, P., Relton, C.L. and Suderman, M.J. Epigenetic modelling of former, current and never smokers. Clinical Epigenetics, 2021, 13, 1-13.
  
  Wen, D., Shi, J., Liu, Y., He, W., Qu, W., Wang, C., Xing, H., Cao, Y., Li, J. and Zha, L. DNA methylation analysis for smoking status prediction in the Chinese population based on the methylation-sensitive single-nucleotide primer extension method. Forensic Science International, 2022, 339, 111412.
  }
  \examples{
    data(melon)
    # note, melon is not a complete dataset, does not work with all methods
    smokp(melon, method="McCartney", sst=NULL)
    
    
  }
