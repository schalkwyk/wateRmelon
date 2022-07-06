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
smokp(betas, method = c('Teschendorff', 'Yu', 'Gao', 'McCartney', 'Yang', 'Sugden', 'Odintsova', 'Zhang', 'Christiansen', 'SSt', 'Packyears', 'Cessation'), ...)
}
\arguments{
  \item{betas}{
Matrix of betas or MethyLumiSet or MethylSet object.
}
  \item{method}{
Currently: 'Teschendorff', 'Yu', 'Gao', 'McCartney', 'Yang', 'Sugden', 'Odintsova', 'Zhang', 'Christiansen', 'SSt', 'Packyears', 'Cessation'. Default is "SSt".
}
  \item{sst}{
Named vector describing smoking status (Current, Former, or Never) of participant for each sample (where names match rownames of betas).
}
%\details{}
\value{
Returns data frame of predicted smoking per sample.
}
\references{
Zhang, Y., Florath, I., Saum, K. U., & Brenner, H. (2016). Self-reported smoking, serum cotinine, and blood DNA methylation. Environmental research, 146, 395-403.
Teschendorff, A.E., Yang, Z., Wong, A., Pipinikas, C.P., Jiao, Y., Jones, A., Anjum, S., Hardy, R., Salvesen, H.B., Thirlwell, C. and Janes, S.M. (2015) Correlation of smoking-associated DNA methylation changes in buccal cells with DNA methylation changes in epithelial cancer. JAMA oncology, 1(4), pp.476-485.
Gao, X., Zhang, Y., Breitling, L.P. and Brenner, H. (2016) Relationship of tobacco smoking and smoking-related DNA methylation with epigenetic age acceleration. Oncotarget, 7(30), p.46878.
Yang, Y., Gao, X., Just, A.C., Colicino, E., Wang, C., Coull, B.A., Hou, L., Zheng, Y., Vokonas, P., Schwartz, J. and Baccarelli, A.A. (2019) Smoking-related DNA methylation is associated with DNA methylation phenotypic age acceleration: The veterans affairs normative aging study. International journal of environmental research and public health, 16(13), p.2356.
Yu, H., Raut, J.R., Schöttker, B., Holleczek, B., Zhang, Y. and Brenner, H. (2020) Individual and joint contributions of genetic and methylation risk scores for enhancing lung cancer risk stratification: data from a population-based cohort in Germany. Clinical epigenetics, 12(1), pp.1-11.
Christiansen, C., Castillo-Fernandez, J.E., Domingo-Relloso, A., Zhao, W., Moustafa, J.E.S., Tsai, P.C., Maddock, J., Haack, K., Cole, S.A., Kardia, S.L.R. and Molokhia, M. (2021) Novel DNA methylation signatures of tobacco smoking with trans-ethnic effects. Clinical epigenetics, 13(1), pp.1-13.
McCartney, D.L., Hillary, R.F., Stevenson, A.J., Ritchie, S.J., Walker, R.M., Zhang, Q., Morris, S.W., Bermingham, M.L., Campbell, A., Murray, A.D. and Whalley, H.C. (2018) Epigenetic prediction of complex traits and death. Genome biology, 19(1), pp.1-11.
Sugden, K., Hannon, E.J., Arseneault, L., Belsky, D.W., Broadbent, J.M., Corcoran, D.L., Hancox, R.J., Houts, R.M., Moffitt, T.E., Poulton, R. and Prinz, J.A. (2019) Establishing a generalized polyepigenetic biomarker for tobacco smoking. Translational psychiatry, 9(1), pp.1-12.
Odintsova, V.V., Rebattu, V., Hagenbeek, F.A., Pool, R., Beck, J.J., Ehli, E.A., van Beijsterveldt, C.E., Ligthart, L., Willemsen, G., De Geus, E.J. and Hottenga, J.J. (2021) Predicting complex traits and exposures from polygenic scores and blood and buccal DNA methylation profiles. Frontiers in Psychiatry, 12.
Elliott, H.R., Tillin, T., McArdle, W.L., Ho, K., Duggirala, A., Frayling, T.M., Smith, G.D., Hughes, A.D., Chaturvedi, N. and Relton, C.L. (2014) Differences in smoking associated DNA methylation patterns in South Asians and Europeans. Clinical epigenetics, 6(1), pp.1-10.
Zhang, Y., Schöttker, B., Florath, I., Stock, C., Butterbach, K., Holleczek, B. & Brenner, H. (2016). Smoking-associated DNA methylation biomarkers and their predictive value for all-cause and cardiovascular mortality. Environmental health perspectives, 124(1), 67-74.
}
\author{
Original Functions: See References.

wateRmelon Implementation: Tyler Gorrie-Stone, Leo Schalkwyk, Louis El Khoury
}
% \note{}

\examples{
data(melon)
smokp(melon, method="all", sst=NULL)
}
