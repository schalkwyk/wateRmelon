# Revision Date: 01-02-2016

outlyx <- function(x, iqr=TRUE, iqrP=2, pc=1,
                   mv=TRUE, mvP=0.15, plot=TRUE, ...) { # {{{
### Computes outliers within methylomic datasets
 # Arguments:
 # x    : Methylumi/Minfi object or raw betas.
 #
 # pc   : The desired principal component for outlier identification
 #
 # iqr  : Logical, to indicate whether to determine outliers by
 #        interquartile ranges.
 #
 # iqrP : The number of interquartile ranges one wishes to
 #        discriminate from the upper and lower quartiles.
 #        Default = 2, Tukey's rule suggests 1.5
 #
 # mv   : Logical, to indicate whether to determine outliers
 #        using distance measures using a modified version of pcout
 #        from mvoutlier.
 #
 # mvP  : The threshold bywhich one wishes to screen
 #        data based on the final weight output from
 #        the pcout function
 #        Value between 0 and 1. Default is 0.15
 # plot : Logical, to indicate if a graphical device to display
 #        sample outlyingness.
 #
###
  # Initialising objects to be used:
  df <- list() # Converted into dataframe later on.
  # Removing probes with NA values.
  x <- na.omit(x)

  if(iqr){ # {{{
    pccompbetx <- prcomp(x, retx=FALSE)
    out1 <- iqrFun(pccompbetx$rot, pc=pc, iqrP=iqrP)
    v1 <- colnames(x) %in% out1[[1]]==TRUE
    df[["iqr"]] <- v1
    low <- min(c(min(pccompbetx$rot[,pc]),out1[["low"]]))-out1[[2]] # Used if Plot=T
    upp <- max(c(max(pccompbetx$rot[,pc]),out1[["hi"]]))+out1[[2]]  # Used if Plot=T
  } # }}}

  if(mv){ # {{{
    tbetx <- t(x)
    out2 <- mvFun(tbetx, mvP=mvP)
    v2 <- colnames(x) %in% out2[[1]]==TRUE
    df[["mv"]] <- v2
  } # }}}

  if(plot&mv&iqr){
    plot(pccompbetx$rot[,pc],
         out2[[2]],
         xlim=c(low, upp),
         xlab="Transformed Betas",
         ylab="Final Weight", ...)
    abline(v=c(out1[["low"]],out1[["hi"]]),
           h=mvP,
           lty=2)
    rect(xleft=low,
         xright=out1[["low"]],
         ybottom=0,
         ytop=mvP,
         col="red",
         density=10)
    rect(xleft=out1[["hi"]],
         xright=upp,
         ybottom=0,
         ytop=mvP,
         col="red",
         density=10)
  }

 # Possibly a better way to do the below.
  df[["outliers"]] <- if(mv){df[["mv"]]==TRUE}&if(iqr){df[["iqr"]]==TRUE}
  df <- data.frame(df)
  rownames(df) <- colnames(x)
  return(df)
} # }}}

iqrFun <- function(x, pc, iqrP){ # {{{
  quantilesbx <- apply(x, 2, quantile) # fivenum also works
  IQR <- quantilesbx[4, pc] - quantilesbx[2, pc]
  thresh <- iqrP*IQR
  hiOutlyx <- quantilesbx[4,pc] + thresh # Upper threshold
  loOutlyx <- quantilesbx[2,pc] - thresh # Lower threshold
  outhi <- x[, pc] > hiOutlyx # Upper Outliers
  outlo <- x[, pc] < loOutlyx # Lower Outliers
  return(list(c(rownames(x)[outhi], rownames(x)[outlo]),
              IQR, "hi" = hiOutlyx, "low" = loOutlyx))
} # }}}

mvFun <- function(x, mvP, ...){ # {{{
  pcoutbetx <- pcouted(x,...)
               #pcout(x)
  outmvbetx <- pcoutbetx < mvP
               #pcoutbetx$wfinal < mvP
  return(list(rownames(as.matrix(pcoutbetx[outmvbetx])),
              pcoutbetx))
} # }}}

pcouted <- function(x,explvar=0.99,crit.M1=1/3,crit.c1=2.5,
                    crit.M2=1/4,crit.c2=0.99,cs=0.25,
                    outbound=0.25, ...){ # {{{
  # Modified version of the pcout function from mvoutlier to
  # calculated distance measures for methylomic data.
  # Two minor things have been changed, mainly how the function
  # copes with x.mad values = 0 and the output.
  #
  #
  x.mad=apply(x,2,mad)
  x <- x[,!x.mad==0] # Cheat to allow it to continue to function.
  x.mad <- x.mad[!x.mad==0]
  p = ncol(x)
  n = nrow(x)
  #
  # PHASE 1:
  # Step 1: robustly sphere the data:
  x.sc <- scale(x,apply(x,2,median),x.mad)
  # Step 2: PC decomposition; compute p*, robustly sphere:
  x.svd <- svd(scale(x.sc,TRUE,FALSE))
  a <- x.svd$d^2/(n-1)
  p1 <- (1:p)[(cumsum(a)/sum(a)>explvar)][1]

  x.pc <- x.sc%*%x.svd$v[,1:p1]
  xpc.sc <- scale(x.pc,apply(x.pc,2,median),apply(x.pc,2,mad))

  # Step 3: compute robust kurtosis weights, transform to distances:
  wp <- abs(apply(xpc.sc^4,2,mean)-3)

  xpcw.sc <- xpc.sc%*%diag(wp/sum(wp))
  xpc.norm <- sqrt(apply(xpcw.sc^2,1,sum))
  x.dist1 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)

  # Step 4: determine weights according to translated biweight:
  M1 <- quantile(x.dist1,crit.M1)
  const1 <- median(x.dist1)+crit.c1*mad(x.dist1)
  w1 <- (1-((x.dist1-M1)/(const1-M1))^2)^2
  w1[x.dist1<M1] <- 1
  w1[x.dist1>const1] <- 0

  #
  # PHASE 2:
  # Step 5: compute Euclidean norms of PCs and their distances:
  xpc.norm <- sqrt(apply(xpc.sc^2,1,sum))
  x.dist2 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)

  # Step 6: determine weight according to translated biweight:
  M2 <- sqrt(qchisq(crit.M2,p1))
  const2 <- sqrt(qchisq(crit.c2,p1))
  w2 <- (1-((x.dist2-M2)/(const2-M2))^2)^2
  w2[x.dist2<M2] <- 1
  w2[x.dist2>const2] <- 0
  #
  # Combine PHASE1 and PHASE 2: compute final weights:
  # Changed output slightly
  wfinal <- (w1+cs)*(w2+cs)/((1+cs)^2)
  return(wfinal)
} # }}}
