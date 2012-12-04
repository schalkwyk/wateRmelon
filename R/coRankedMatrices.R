coRankedMatrices <-
function(matrix1, matrix2){

	rownames1 <- rownames(matrix1)
	rownames2 <- rownames(matrix2)

	matrix1 <- matrix1[sort(rownames1, index.return=TRUE)$ix,]
	matrix2 <- matrix2[sort(rownames2, index.return=TRUE)$ix,]
	
	equalNames <- rownames(matrix1) == rownames(matrix2)
	indexFalse <- which(equalNames==FALSE)
	if(length(indexFalse)>0){
		print("WARNING: matrices row names do not have same rank or matrices row names differ. Check the 'return' result for mismatchs index...")
		return(indexFalse)
	}

	return(list(matrix1=matrix1, matrix2=matrix2))
}
