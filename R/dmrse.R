#' Standard error of iDMR 450k array DNA methylation features
#' 
#' Imprinting differentially methylated regions (iDMRs) are expected to be
#' approximately half methylated, as is observed at the 227 probes in known
#' iDMRs. These functions calculate measures of dispersion for the beta values
#' at these CpG sites, of which the most useful is dmrse_row, which is the
#' between-sample standard error.
#' 
#' 
#' @aliases dmrse dmrse_col dmrse_row
#' @param betas a matrix of betas (default method), a \code{MethyLumiSet}
#' object (\code{methylumi} package), a \code{MethylSet} or \code{RGChannelSet}
#' object ( \code{minfi} package) or a \code{exprmethy450} object (\code{IMA}
#' package).
#' @param idmr %% ~~Describe \code{idmr} here~~ a character vector of iDMR
#' probe names such as returned by iDMR()
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ... return a standard error of the mean of betas for all samples
#' and iDMR probes (dmrse) or the standard error of the mean for just the
#' between sample component(dmrse_row) or between probe(dmrse_col) component.
#' @author Leo Schalkwyk <lschal@@essex.ac.uk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link{seabi}}, a sex-difference metric, and \code{\link{genki}}, based
#' on SNPs.
#' @references Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC: A
#' data-driven approach to preprocessing Illumina 450K methylation array data
#' (submitted)
#' @examples
#' 
#' 
#'   #MethyLumiSet method
#'      data(melon)
#'      dmrse(melon)
#' 
#'   #MethyLumiSet method after normalization
#'      melon.dasen <- dasen(melon)
#'      dmrse(melon.dasen)
#' 
#' @export dmrse
dmrse <-
function(betas, idmr=iDMR()) {  # formerly SDO - both
    idmr <- idmr[idmr %in% rownames(betas)]
    sd(as.numeric(betas[idmr,]), na.rm=T) / sqrt(dim(betas)[2])
}


iDMR <- function (){


c(
"cg00000029", "cg00155882", "cg00576435", "cg00702231", "cg00765653", 
"cg00766368", "cg00906934", "cg01026744", "cg01071811", "cg01174175", 
"cg01355739", "cg01466133", "cg01578057", "cg01585333", "cg01784351", 
"cg01873334", "cg01893176", "cg02146091", "cg02162069", "cg02219360", 
"cg02472486", "cg02611863", "cg02722214", "cg02864690", "cg03061677", 
"cg03085377", "cg03384175", "cg03401726", "cg03422070", "cg03562868", 
"cg03615235", "cg03803002", "cg04489586", "cg04521244", "cg04810287", 
"cg04975775", "cg04995311", "cg05065100", "cg05096321", "cg05326984", 
"cg05556276", "cg05558390", "cg05740879", "cg05816130", "cg05862114", 
"cg05898246", "cg06000530", "cg06163629", "cg06179486", "cg06288089", 
"cg06617468", "cg06695761", "cg06982169", "cg07156273", "cg07189342", 
"cg07360805", "cg07539802", "cg07595203", "cg07643061", "cg07878171", 
"cg08263357", "cg08338216", "cg08359167", "cg08402058", "cg08446215", 
"cg08633479", "cg08835688", "cg08962841", "cg09240436", "cg09452478", 
"cg09512080", "cg09518720", "cg09541000", "cg09701145", "cg10140536", 
"cg10154633", "cg10204755", "cg10642330", "cg10981598", "cg11174847", 
"cg11175683", "cg11297256", "cg11300971", "cg11399589", "cg11532302", 
"cg11562309", "cg11613559", "cg11666921", "cg11948874", "cg11985632", 
"cg11993252", "cg12077660", "cg12205903", "cg12243267", "cg12298755", 
"cg12532169", "cg12757684", "cg12759554", "cg12862537", "cg12903171", 
"cg13431205", "cg13591710", "cg13610072", "cg13690564", "cg13708635", 
"cg13713522", "cg13790727", "cg13797824", "cg13828758", "cg13901453", 
"cg13960339", "cg14161241", "cg14175568", "cg14203179", "cg14243741", 
"cg14306330", "cg14392746", "cg14469070", "cg14597908", "cg14728235", 
"cg14765818", "cg14849423", "cg14873490", "cg14958441", "cg15302378", 
"cg15330298", "cg15473473", "cg15651941", "cg15777825", "cg15815607", 
"cg15886040", "cg16303279", "cg16326421", "cg16463460", "cg16492735", 
"cg16547341", "cg16574793", "cg16648571", "cg16739686", "cg17580798", 
"cg17635114", "cg17643025", "cg17769238", "cg17840843", "cg17895149", 
"cg18083595", "cg18433380", "cg18481241", "cg18506672", "cg18602919", 
"cg18607468", "cg18619398", "cg18723276", "cg18898992", "cg19151808", 
"cg19177307", "cg19296354", "cg19344806", "cg19427472", "cg19617948", 
"cg19663555", "cg19782411", "cg20252111", "cg20365618", "cg20479660", 
"cg20601276", "cg20699737", "cg20783699", "cg21137515", "cg21200654", 
"cg21526238", "cg21588305", "cg21625881", "cg21771834", "cg21926091", 
"cg21952820", "cg22259242", "cg22298088", "cg22378853", "cg22421148", 
"cg22502459", "cg22510412", "cg22551578", "cg22694067", "cg22807877", 
"cg22849953", "cg22943498", "cg23342787", "cg23401210", "cg23460430", 
"cg23566503", "cg23605670", "cg23714917", "cg23742152", "cg23757721", 
"cg23954636", "cg24065722", "cg24338351", "cg24409677", "cg24564503", 
"cg24675557", "cg24762053", "cg24975842", "cg24977055", "cg25306939", 
"cg25437674", "cg25645178", "cg25712981", "cg25962605", "cg26083330", 
"cg26094482", "cg26104781", "cg26349266", "cg26406256", "cg26469586", 
"cg26503018", "cg26547719", "cg26725891", "cg26908876", "cg27001184", 
"cg27120649", "cg27150681", "cg27216384", "cg27323091", "cg27372170", 
"cg27409910", "cg27589003"
)

}
