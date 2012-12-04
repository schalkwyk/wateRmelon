getMethylumiBeta <-
function(methylumi){

	fNames <- featureNames(methylumi)
	u <- unmethylated(methylumi)
	m <- methylated(methylumi)

	rm(methylumi)

	#check and "correct" for negative values
	indexNegU <- which(u < 0, arr.ind=TRUE)
	indexNegM <- which(m < 0, arr.ind=TRUE)
	u[indexNegU] <- 0
	m[indexNegM] <- 0

	rm(indexNegU, indexNegM)

	beta <- m/(m + u + 100)
	rm(u, m)

	rownames(beta) <- fNames

	return(beta)
}
