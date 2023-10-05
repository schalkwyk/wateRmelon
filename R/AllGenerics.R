# generic cornucopia  S4

setGeneric ( name= "naten"          )
setGeneric ( name= "betaqn"         )
setGeneric ( name= "nanet"          )
setGeneric ( name= "nanes"          )
setGeneric ( name= "nasen"          )
setGeneric ( name= "danes"          )
setGeneric ( name= "danet"          )
setGeneric ( name= "daten1"         )
setGeneric ( name= "daten2"         )
setGeneric ( name= "dasen"          )
setGeneric ( name= "danen"          )
setGeneric ( name= "tost"           )
setGeneric ( name= "fuks"           )
setGeneric ( name= "swan"           )
setGeneric ( name= "colnames"       )
setGeneric ( name= "pfilter"        )
setGeneric ( name= "genki"          )
setGeneric ( name= "seabi"          )
setGeneric ( name= "dmrse"          )
setGeneric ( name= "dmrse_row"      )
setGeneric ( name= "dmrse_col"      )
setGeneric ( name= "BMIQ"           )
setGeneric ( name= "bscon"          )
#setGeneric ( name= "as.methylumi"   )
setGeneric ( name= "outlyx"         )
setGeneric ( name= "pwod"           )
setGeneric ( name= "agep"           )
setGeneric ( name= "uSexQN"         )
setGeneric ( name= "estimateCellCounts.wmln" )


# S3 generics to work with gds (bigmelon)

epicv2clean <- function(x){
   UseMethod('epicv2clean', x)
}
