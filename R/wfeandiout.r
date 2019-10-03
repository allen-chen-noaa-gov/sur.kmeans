wfeandiout <- function(FullTableFull,pooledvcov,testvars,nontestvars,SURbetaFE,
    YYlist,XXlist) {
    #' WFE estimates for construction of dispersion test
    #'
    #' WFE estimates for use in constructing dispersion test, gets called as
    #' part of constructing the test in kmeans_inner
    #'
    #' @param FullTableFull data
    #' @return Weighted fixed effects estimates (and associated sigma weights)
    #' @export
    #' @examples
    #'
    
    FullTableFull$Group <- as.factor(FullTableFull$Group)
    FullTableFull$Vessel <- as.factor(FullTableFull$Vessel)
    FullTableFull$EQ <- as.numeric(FullTableFull$EQ)
    
    betaWFE <- list()
    betai <- list()
    sigmai <- list()
    residsend <- list()
    
    for (k in levels(FullTableFull$Group)) {
        
        FullTable2 <- subset(FullTableFull, Group == k)
        FullTable2$Vessel <- droplevels(FullTable2$Vessel)
        
        WFEpart1 <- list()
        WFEpart2 <- list()
        
        for (i in levels(FullTable2$Vessel)) {
            
            FullTable <- subset(FullTable2, Vessel == i)
            
            XX <- as.matrix(FullTable[unique(unlist(XXlist))])
            
            YY <- as.matrix(rowSums(FullTable[unique(unlist(YYlist))]))
            
            pooledBigSigma <- kronecker(solve(pooledvcov[[k]]),
			    diag((dim(YY)[1]/4)))
            cholPBS <- chol(pooledBigSigma)
            
			#note here that t(cholPBS %*% XX) %*% (cholPBS %*% XX) and
			    #t(XX) %*% (pooledBigSigma) %*% (XX) are numerically identical
				#(e.g. try all.equal()).
				#If applying FWL the order matters however, and need to
				#premultiply first.
            PcholXX <- cholPBS %*% XX
            PcholYY <- cholPBS %*% YY
            
			#see remarks 3 and 4 for form of partitioned unbalanced dispersion 
			    #tests in Pesaran & Yamagata 2008
            PcholXXtest <- PcholXX[, testvars]
            PcholXXnontest <- PcholXX[, nontestvars]
            PcholPartM <- diag(dim(YY)[1]) - PcholXXnontest %*%
			    solve(t(PcholXXnontest) %*% PcholXXnontest) %*%
				t(PcholXXnontest)
            
            SURbetaFEtest <- SURbetaFE[[k]]
            
            residsend[[k]] <- t(PcholYY - PcholXXtest %*% SURbetaFEtest) %*%
			    (PcholYY - PcholXXtest %*% SURbetaFEtest)
            
            tempsigma <- (t(PcholPartM %*% PcholYY - PcholPartM %*%
			    PcholXXtest %*% SURbetaFEtest) %*%(PcholPartM %*% PcholYY -
				PcholPartM %*% PcholXXtest %*% SURbetaFEtest))/
				(dim(YY)[1] - (length(nontestvars)) - 1)
            
            WFEpart1[[i]] <- (t(PcholPartM %*% PcholXXtest) %*% (PcholPartM %*%
                PcholXXtest))/as.numeric(tempsigma)
            WFEpart2[[i]] <- (t(PcholPartM %*% PcholXXtest) %*% (PcholPartM %*% 
                PcholYY))/as.numeric(tempsigma)
            
            yvars <- list()
            ones <- list()
            xvars <- list()
            singlebeta <- list()
            resideq <- list()
            
            for (a in 1:max(FullTable$EQ)) {
                
                xtemp <- FullTable[FullTable$EQ == a, ]
                
                yvars[[a]] <- as.matrix(xtemp[YYlist[a]])
                xvars[[a]] <- xtemp[XXlist[[a]]]
                
                ones[[a]] <- diag(dim(yvars[[a]])[1]) -
				    rep(1, dim(xtemp)[1]) %*%
                    solve(t(rep(1, dim(xtemp)[1])) %*%
					rep(1, dim(xtemp)[1])) %*%
                    t(rep(1, dim(xtemp)[1]))
                
                xvars[[a]] <- xvars[[a]][, colSums(xvars[[a]] == 1) !=
				    dim(xvars[[a]])[1]]
                xvars[[a]] <- as.matrix(xvars[[a]][, colSums(xvars[[a]] == 0) !=
                    dim(xvars[[a]])[1]])
                
                singlebeta[[a]] <- solve(t(xvars[[a]]) %*% ones[[a]] %*%
				    xvars[[a]]) %*% (t(xvars[[a]]) %*% ones[[a]] %*% yvars[[a]])
                
                resideq[[a]] <- ones[[a]] %*% yvars[[a]] - ones[[a]] %*%
				    xvars[[a]] %*% singlebeta[[a]]
                
                resideq[[a]] <- resideq[[a]] - mean(resideq[[a]])
                
            }
            
            C <- ((dim(FullTable)[1]/max(FullTable$EQ)) - 1)^(-1) *
			    t(matrix(unlist(resideq), dim(FullTable)[1]/max(FullTable$EQ),
				max(FullTable$EQ))) %*% matrix(unlist(resideq),
				dim(FullTable)[1]/max(FullTable$EQ), max(FullTable$EQ))
            
            singlevcov <- C
            
            singleBigSigma <- kronecker(solve(singlevcov), diag((dim(YY)[1]/4)))
            cholSBS <- chol(singleBigSigma)
            
            cholXX <- as.data.frame(cholSBS %*% XX)
            cholXXtest <- as.matrix(cholXX[, testvars])
            cholXXnontest <- as.matrix(cholXX[, nontestvars])
            
            cholPartM <- diag(dim(YY)[1]) - cholXXnontest %*%
			    solve(t(cholXXnontest) %*% cholXXnontest) %*% t(cholXXnontest)
            cholYY <- cholSBS %*% YY
            
            singlebeta <- solve(t(cholPartM %*% cholXXtest) %*%
			    (cholPartM %*% cholXXtest)) %*% (t(cholPartM %*% cholXXtest) %*%
				(cholPartM %*% cholYY))
            
            sigmai[[i]] <- tempsigma
            betai[[i]] <- singlebeta
            
        }
        
        betaWFE[[k]] <- solve(Reduce("+", WFEpart1)) %*% Reduce("+", WFEpart2)
        
    }
    
    return(list(betai = betai, betaWFE = betaWFE, sigmai = sigmai,
	    residsend = residsend))
    
}
