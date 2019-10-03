pooledvcovout <- function(FullTableFull,betaFE,YYlist,XXlist) {
    #' Estimate variance-covariance matrix
    #'
    #' Estimate variance-covariance matrix
    #'
    #' @param FullTableFull data
    #' @return Variance-covariance matrix
    #' @export
    #' @examples
    #'
    
    FullTableFull$Group <- as.factor(FullTableFull$Group)
    FullTableFull$Vessel <- as.factor(FullTableFull$Vessel)
    FullTableFull$EQ <- as.numeric(FullTableFull$EQ)
    
    pooledvcov <- list()
    
    for (k in levels(FullTableFull$Group)) {
        
        FullTable2 <- subset(FullTableFull, Group == k)
        FullTable2$Vessel <- droplevels(FullTable2$Vessel)
        FullTable <- FullTable2
        
        yvars <- list()
        ones <- list()
        xvars <- list()
        FEpart1 <- list()
        FEpart2 <- list()
        dummymatrix <- list()
        resideq <- list()
        
        for (a in 1:max(FullTable$EQ)) {
            
            xtemp <- FullTable[FullTable$EQ == a, ]
            
            yvars[[a]] <- as.matrix(xtemp[YYlist[a]])
            xvars[[a]] <- xtemp[XXlist[[a]]]
            
            ones[[a]] <- diag(dim(yvars[[a]])[1]) - rep(1, dim(xtemp)[1]) %*%
			    solve(t(rep(1, dim(xtemp)[1])) %*% rep(1, dim(xtemp)[1])) %*%
				t(rep(1, dim(xtemp)[1]))
            
            xvars[[a]] <- xvars[[a]][, colSums(xvars[[a]] == 1) !=
			    dim(xvars[[a]])[1]]
            xvars[[a]] <- as.matrix(xvars[[a]][, colSums(xvars[[a]] == 0) !=
			    dim(xvars[[a]])[1]])
            
            FEpart1[[a]] <- t(xvars[[a]]) %*% ones[[a]] %*% xvars[[a]]
            FEpart2[[a]] <- t(xvars[[a]]) %*% ones[[a]] %*% yvars[[a]]
            
            dummymatrix <- tryCatch(model.matrix(~Vessel - 1, data = xtemp),
			    error = function(e) rep(1, dim(yvars[[a]])[1]))
            BigM <- diag(dim(yvars[[a]])[1]) - dummymatrix %*%
			    solve(t(dummymatrix) %*% dummymatrix) %*% t(dummymatrix)
            
            resideq[[a]] <- BigM %*% yvars[[a]] - BigM %*% xvars[[a]] %*%
			    betaFE[[k]][[a]]
            
            resideq[[a]] <- resideq[[a]] - mean(resideq[[a]])
            
        }
        
        C <- ((dim(FullTable)[1]/max(FullTable$EQ)) - 1)^(-1) *
		    t(matrix(unlist(resideq), dim(FullTable)[1]/max(FullTable$EQ),
			max(FullTable$EQ))) %*% matrix(unlist(resideq),
            dim(FullTable)[1]/max(FullTable$EQ), max(FullTable$EQ))
        
        pooledvcov[[k]] <- C
        
    }
    
    return(pooledvcov)
    
}
