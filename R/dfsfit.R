dfsfit <-
function(
   mn, 
   onetwo,
   roco = substring(colnames(mn), regexpr("R0[1-9]C0[1-9]", colnames(mn))), # Robust solution
   ...
){
   mdf<-apply(mn,2,dfs2,onetwo)

   if (! is.null(roco) ) {
      scol  <- as.numeric(substr(roco,6,6))
      srow  <- as.numeric(substr(roco,3,3))
      fit   <- try(  lm( mdf ~ srow + scol ), silent=TRUE) 
      if (! inherits (fit, "try-error") ) {mdf   <- fit$fitted.values}
      else { message ('Sentrix position model failed, skipping') }
   }
   otcor <-  matrix(
      rep( mdf, sum(onetwo=='I')),
      byrow=T, 
      nrow=sum(onetwo=='I')
   )
   mn[onetwo=='I',] <- mn[onetwo=='I',] - otcor
   mn
}
