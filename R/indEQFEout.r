indEQFEout <- function(FullTableFull,YYlist,XXlist) {
    #' Individual FE estimates for construction of dispersion test
    #'
    #' Individual FE estimates for use constructing in dispersion test, gets 
    #' called as part of constructing the test in kmeans_inner
	#'
    #' @param FullTableFull data
    #' @return Parameter estimates
    #' @export
    #' @examples
    #'
    
    FullTableFull$Group <- as.factor(FullTableFull$Group)
    FullTableFull$Vessel <- as.factor(FullTableFull$Vessel)
    FullTableFull$EQ <- as.numeric(FullTableFull$EQ)
    
    betaFE <- list()
    
    for (k in levels(FullTableFull$Group)) {
        
        FullTable2 <- subset(FullTableFull, Group == k)
        FullTable2$Vessel <- droplevels(FullTable2$Vessel)
        
        FEpart1G <- list()
        FEpart2G <- list()
        
        for (i in levels(FullTable2$Vessel)) {
            
            FullTable <- subset(FullTable2, Vessel == i)
            
            yvars <- list()
            ones <- list()
            xvars <- list()
            FEpart1 <- list()
            FEpart2 <- list()
            
            for (a in 1:max(FullTable$EQ)) {
                
                xtemp <- FullTable[FullTable$EQ == a, ]
                
                yvars[[a]] <- as.matrix(xtemp[YYlist[a]])
                xvars[[a]] <- xtemp[XXlist[[a]]]
                
                ones[[a]] <- diag(dim(yvars[[a]])[1]) -
				    rep(1, dim(xtemp)[1]) %*% solve(t(rep(1, dim(xtemp)[1])) %*%
                    rep(1, dim(xtemp)[1])) %*% t(rep(1, dim(xtemp)[1]))
                
                xvars[[a]] <- xvars[[a]][, colSums(xvars[[a]] == 1) !=
				    dim(xvars[[a]])[1]]
                xvars[[a]] <- as.matrix(xvars[[a]][, colSums(xvars[[a]] == 0) !=
                    dim(xvars[[a]])[1]])
                
                FEpart1[[a]] <- t(xvars[[a]]) %*% ones[[a]] %*% xvars[[a]]
                FEpart2[[a]] <- t(xvars[[a]]) %*% ones[[a]] %*% yvars[[a]]
                
            }
            
            FEpart1G[[i]] <- FEpart1
            FEpart2G[[i]] <- FEpart2
            
        }
        
        betaFEG <- list()
        
        for (a in 1:max(FullTable$EQ)) {
            
            betaFEG[[a]] <- solve(Reduce("+", lapply(FEpart1G, "[[", a))) %*%
			    Reduce("+", lapply(FEpart2G, "[[", a))
            
        }
        
        betaFE[[k]] <- betaFEG
        
    }
    
    return(betaFE)
    
}
