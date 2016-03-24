### Bisulfite conversion QC function, derived from E. Hannon's
### script to calculate bisulfite conversion values
### Louis el Khoury and Leo Schalkwyk 2015

# see also bscon_methy and bscon_minfi

bscon <- function(x, ...) { UseMethod (bscon, x ) }



