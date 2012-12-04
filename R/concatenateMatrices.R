concatenateMatrices <-
function(matrix1, matrix2){
		
		commonRows <- intersect(rownames(matrix1), rownames(matrix2))
		if(length(commonRows) > 0){
			matrix1 <- matrix1[which(is.element(rownames(matrix1), commonRows)),]
			matrix2 <- matrix2[which(is.element(rownames(matrix2), commonRows)),]
		
			matrix1 <- matrix1[sort(rownames(matrix1), index.return=TRUE)$ix,]
			matrix2 <- matrix2[sort(rownames(matrix2), index.return=TRUE)$ix,]
		
			return(cbind(matrix1, matrix2))
		}
		else{
			print("WARNING ! in 'concatenateMatrices', matrix1 and matrix2 must share at least a common row name.")
			return("WARNING ! in 'concatenateMatrices', matrix1 and matrix2 must share at least a common row name")
		}
}
