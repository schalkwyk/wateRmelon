
# methods for  exprmethy450  (IMA)

#setClass("exprmethy450",
# representation("VIRTUAL",
# description = "character"))


# betaqn <- function (bn){ 
setMethod(
   f= "betaqn",
   signature(bn="exprmethy450"),
   definition=function(bn){
   if(!library(IMA, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load IMA package')
      }
         object <- bn
         object@bmatrix <- betaqn (
         bn = object@bmatrix 
         ) 
      object
   }
)
## exprmethy450 only has betas ##
# naten <- function (mn, un, fudge=100 ) {
#nanet <- function(mn, un, fudge=100){
#nanes <- function(mns, uns, onetwo, fudge=100){
#danes <- function(mn, un, onetwo, fudge=100, ...){ 
#danet <- function(mn, un, onetwo, fudge=100, ...){
#daten1 <- function(mn, un, onetwo, fudge=100, ...){
#daten2 <- function(mn, un, onetwo, fudge=100, ...){
#ot <- function ( mns, uns, onetwo, fudge=100) {
#dasen <- function(mns, uns, onetwo, fudge=100, ...){
#danen <- function ( mns, uns, onetwo, fudge=100, ...){
#tost <- function( mn, un, da, pn ) {
#peak.correction <- function (data, anno) {
setMethod(
   f= "fuks",
   signature(data="exprmethy450"),
   definition=function(data, anno=NULL){
      if(!library(IMA, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load IMA package')
      }
      object <- data
      object@bmatrix <- fuks (
         data = object@bmatrix, 
         anno = object@anno
      ) 
      object
   }
)

#swan <- function (mn, un, qc ) {

#genki <- function(bn, g=getsnp(rownames(bn)), se=TRUE ){
setMethod(
   f= "genki",
   signature(bn="exprmethy450"),
   definition=function(bn, g=NULL, se=TRUE){
      if(!library(IMA, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load IMA package')
      }
      object <- bn
      bn     <- object@bmatrix
      g      <- getsnp(rownames(bn))
      genki( bn, g, se ) 
     
   }
)

#dmrse <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse",
   signature(betas="exprmethy450"),
   definition=function(betas, idmr=iDMR()){
      if(!library(IMA, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load IMA package')
      }
      object <- betas
      betas     <- object@bmatrix
      dmrse( betas, idmr ) 
     
   }
)


#dmrse_row <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse_row",
   signature(betas="exprmethy450"),
   definition=function(betas, idmr=iDMR()){
      if(!library(IMA, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load IMA package')
      }
      object <- betas
      betas     <- object@bmatrix
      dmrse_row( betas, idmr ) 
     
   }
)

#dmrse_col <- function(betas, idmr=iDMR)
setMethod(
   f= "dmrse_col",
   signature(betas="exprmethy450"),
   definition=function(betas, idmr=iDMR()){
      if(!library(IMA, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load IMA package')
      }
      object <- betas
      betas     <- object@bmatrix
      dmrse_col( betas, idmr ) 
     
   }
)

#seabi <- function (bn, stop, sex, X){
setMethod(
   f= "seabi",
   signature(bn="exprmethy450"),
   definition=function( bn, stop=1, sex, X ){
      if(!library(IMA, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load IMA package')
      }
      bn     <- object@bmatrix
      seabi( bn, stop, sex, X ) 
     
   }
)

