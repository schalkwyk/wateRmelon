getsnp <-
function(nam) { # takes a vector of rownames

   g <- grep('^rs',nam)
   print ( message ( length(g) , 'SNP data rows found\n' ))
   g
}
