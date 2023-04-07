

admissible_subsets_CV <- function(Y, X, N = 100, steps = 2e2, burnin = 1e2, epsilonGrid, nFold = 5, 
                                  Efficient=FALSE, useMCsample4h=FALSE, addoffset=FALSE, parallel=TRUE, n_cores=4, 
                                  trueModel, MLasso=TRUE){
  
  
  require(foreach); require(MASS);
  
  if (MLasso){
    MLasso_Soln     <- InitialValue_MLasso(Y.m=t(Y), X.m=t(X))
    PropWeights     <- unname(rowSums(MLasso_Soln^2))
  }else{
    MGLasso_Soln    <- quiet(InitialValue_MGLasso(Y.m=t(Y), X.m=t(X), lamG=seq(0.05,4,0.1), fold=10))
    PropWeights     <- rowSums(MGLasso_Soln^2)
  }
  
  
  nnzweights      <- which(PropWeights != 0)
  
  print("MSG Lasso Chosen Model")
  print(nnzweights)
  
  n <- ncol(Y); p <- nrow(X); q <- nrow(Y); 
  
  perm.n          <- sample(1:n, n)
  Y               <- Y[, perm.n]
  X               <- X[, perm.n]
  
  fold            <- cut(seq(1,n),breaks=nFold,labels=FALSE)
  fold_freq       <- table(fold)
  folds.occ       <- as.numeric(names(fold_freq))
  fold_prop       <- as.vector(fold_freq)/n
  
  nFold           <- length(folds.occ)
  
  print(paste0("Number of folds: ", nFold))
  
  print("Fold Proportions")
  print(fold_prop)
  
  pred.MSE.OLS    <- matrix(NA, length(epsilonGrid), nFold)
  pred.MSE.BMA    <- matrix(NA, length(epsilonGrid), nFold)
  
  for (k in 1:nFold){
    train         <- fold != folds.occ[k]
    test          <- fold == folds.occ[k]
    
    Y.train       <- Y[,train]
    X.train       <- X[,train]
    
    print(paste0("Fold: ", k))
    print(paste0("FoldValue: ", folds.occ[k]))
    
    if (parallel){
      
      library(doParallel)
      library(doSNOW)
      cl                  <- makeSOCKcluster(n_cores)
      registerDoSNOW(cl)
      progress            <- function(nfin, tag) { cat(sprintf('tasks completed: %d; tag: %d\n', nfin, tag)) }
      opts                <- list(progress=progress)
      packages_req        <- c("matrixsampling", "tidyverse", "Rfast", "MSGLasso", "MASS", "glmnet")
      
      if (Efficient){
        result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
                                     .export = ls(globalenv()) ) %dopar% 
                            { admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
                                                           PropWeights=PropWeights, epsilon=epsi, true.Sigma = NULL, 
                                                           addoffset=addoffset, useMCsample4h=useMCsample4h,
                                                           trueModel = trueModel, MLasso=MLasso)  }
      } else{
        result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
                                     .export = ls(globalenv()) ) %dopar% 
                                    { admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
                                                         PropWeights=PropWeights, epsilon=epsi, true.Sigma = NULL, 
                                                         addoffset=addoffset, useMCsample4h=useMCsample4h,
                                                         trueModel = trueModel, MLasso=MLasso)  }
      }
      
      parallel::stopCluster(cl)
      
    } else{
      
      if (Efficient){
        result <- list()
        for (i in seq_along(epsilonGrid)){
          result[[i]] <- admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
                                                      PropWeights=PropWeights, epsilon=epsilonGrid[i], true.Sigma = NULL, 
                                                      addoffset=addoffset, useMCsample4h=useMCsample4h, 
                                                      trueModel = trueModel, MLasso=MLasso)
        }
        
      } else{
        result <- list()
        for (i in seq_along(epsilonGrid)){
          result[[i]] <- admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
                                            PropWeights=PropWeights, epsilon=epsilonGrid[i], true.Sigma = NULL, 
                                            addoffset=addoffset, useMCsample4h=useMCsample4h, 
                                            trueModel = trueModel, MLasso=MLasso)
        }
        
      }
    }
    
    print(paste0("Number of test samples in fold: ", k))
    print(sum(test))
    
    for (i in seq_along(epsilonGrid)){
      
      Opt.model              <- result[[i]]$Opt.model
      Opt.B                  <- result[[i]]$Opt.B
      Avg.B                  <- result[[i]]$Avg.B
      Avg.size               <- result[[i]]$Avg.size
      
      print(paste0("For epsilon = ", epsilonGrid[i]))
      cat("Chosen model:", Opt.model)
      cat("\n")
      
      pred.OLS             <- Opt.B %*% X[, test, drop=FALSE] 
      pred.BMA             <- Avg.B %*% X[, test, drop=FALSE]
      
      pred.MSE.BMA[i,k]    <- mean((Y[, test, drop=FALSE] - pred.BMA)^2)
      pred.MSE.OLS[i,k]    <- mean((Y[, test, drop=FALSE] - pred.OLS)^2)
      
    }
    
  }
  
  print("Fold calculations at a glance")
  # print(pred.MSE.BMA); 
  print(pred.MSE.OLS);
  
  MSE.BMA             <- apply(pred.MSE.BMA, 1, weighted.mean, w=fold_prop, na.rm = FALSE)
  MSE.OLS             <- apply(pred.MSE.OLS, 1, weighted.mean, w=fold_prop, na.rm = FALSE)
  
  MSE.BMA.eps         <- epsilonGrid[which.min(MSE.BMA)]
  MSE.OLS.eps         <- epsilonGrid[which.min(MSE.OLS)]

  
  optimal.eps         <- c(MSE.OLS.eps, MSE.BMA.eps)
  names(optimal.eps)  <- c("MSE.OLS", "MSE.BMA")
  
  print(MSE.OLS);
  
  list("optimal.eps" =  optimal.eps, "MSE"=cbind(MSE.OLS, MSE.BMA))
  
}
