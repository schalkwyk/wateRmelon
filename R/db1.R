#' Internal wateRmelon functions for calculating betas
#' 
#' db1 is used for quantile normalizing methylated together with unmethylated
#' (dye bias methods nanet, nanes, danes and danet.  dfs* functions are used
#' for smoothing the background equalization in methods whose names start with
#' d (daten etc).
#' 
#' db1 - quantile normalizes methylated against unmethylated (basic function
#' for dyebuy* dye bias methods). dfsfit - corrects the difference in
#' backgrounds between type I and type II assays and fits a linear model to
#' Sentrix rows and columns if these are available to improve precision where
#' there is a background gradient. dfs2 - finds the difference between type I
#' and type II assay backgrounds for one or more samples. %% ~~ If necessary,
#' more details than the description above ~~
#' 
#' @aliases db1 dfs2 dfsfit
#' @param mn,x matrix of methylated signal intensities, each column
#' representing a sample (default method), or an object for which a method is
#' available. For dfsfit and dfs2 this can also be a matrix of unmethylated
#' intensities.
#' @param un matrix of unmethylated signal intensities, each column
#' representing a sample (default method) or NULL when mn is an object
#' containing methylated and unmethylated values.
#' @param onetwo character vector or factor of length nrow(mn) indicating assay
#' type 'I' or 'II'
#' @param roco roco for dfsfit giving Sentrix rows and columns.  This allows a
#' background gradient model to be fit.  This is split from data column names
#' by default.  roco=NULL disables model fitting (and speeds up processing),
#' otherwise roco can be supplied as a character vector of strings like
#' 'R01C01' (3rd and 6th characters used).
#' @return db1 - a list of 2 matrices of intensities, methylated and
#' unmethylated dfsfit - a matrix of adjusted intensities dfs2 - a background
#' offset value
#' 
#' %% ~Describe the value returned %% If it is a LIST, use %% \item{comp1
#' }{Description of 'comp1'} %% \item{comp2 }{Description of 'comp2'} %% ...
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>
#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted)
#' @export db1
db1 <-
function(mn, un ){
   stopifnot(dim(un) == dim(mn))
   a <- dim(un)[2]
   mun <- normalizeQuantiles(cbind(mn,un))
   mnn <- mun[,1:a]
   unn <- mun[,(a+1):(2*a)]
   list(mnn,unn)
}
