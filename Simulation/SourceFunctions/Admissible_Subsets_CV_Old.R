
# In the EAS methodology, epsilon acts a tuning parameter. This function 
# is written to chose the optimal value of epsilon through k-fold 
# cross-validation with respect to RMSE so that we can call Admissible_Subsets.R
# function to get the final model with optimally chosen epsilon

# admissible_subsets_CV <- function(Y, X, N = 100, steps = 2e2, burnin = 1e2, epsilon, nFold = 5, Efficient=FALSE, useMCsample4h=FALSE){
#   
#   n <- ncol(Y); p <- nrow(X); q <- nrow(Y); 
#   fold         <- sample(1:nFold, n, replace = TRUE)
#   Y.pred.GFI   <- matrix(NA,q,n)
#   Y.pred.OLS   <- matrix(NA,q,n)
#   pred.BIC.OLS.fold <- double(nFold)
#   pred.BIC.GFI.fold <- double(nFold)
#   
#   for (k in 1:nFold){
#     train         <- fold != k
#     test          <- fold == k
#     if (Efficient){
#       result        <- admissible_subsets_Efficient(Y[,train],X[,train],N=N,steps=steps,burnin=burnin,PropWeights=NULL,
#                                                     epsilon=epsilon, true.Sigma = NULL, addoffset=TRUE, useMCsample4h=useMCsample4h)
#     } else{
#       result        <- admissible_subsets(Y[,train],X[,train],N=N,steps=steps,burnin=burnin,PropWeights=NULL,
#                                                     epsilon=epsilon, true.Sigma = NULL, addoffset=TRUE, useMCsample4h=useMCsample4h)
#     }
#     
#     models        <- lapply(1:(steps-burnin), function(x) which(result$chain[x,]==1))
#     models.char   <- unlist(lapply(models, paste, collapse = " "))
#     m             <- models[[which(models.char == names(which.max(table(models.char))))[1]]]
#     
#     print(paste0(rep("*",100), collapse = ""))
#     print(paste0("For epsilon = ", epsilon, " and fold = ",k, ", the frequency distribution of the selected models"))
#     print(table(models.char))
#     print(paste0(rep("*",100), collapse = ""))
#     print(paste("Chosen Model at fold", k))
#     print(m)
#     
#     pred.GFI            <- Predict_MLM(Y.train=Y[, train], X.train=X[m, train, drop=FALSE], X.test=X[m, test, drop=FALSE], method="GFI", nsamp=1e3)$yhat
#     pred.OLS            <- Predict_MLM(Y.train=Y[, train], X.train=X[m, train, drop=FALSE], X.test=X[m, test, drop=FALSE], method="OLS", nsamp=1e3)$yhat
#     
#     Y.pred.GFI[,test]   <- pred.GFI
#     Y.pred.OLS[,test]   <- pred.OLS
#     
#     pred.MSE.GFI.fold   <- mean((Y[, test] - pred.GFI)^2)
#     pred.MSE.OLS.fold   <- mean((Y[, test] - pred.OLS)^2)
#     
#     pred.BIC.GFI.fold[k]   <- sum(test)*log(pred.MSE.GFI.fold) + length(m)*q*log(sum(test))
#     pred.BIC.OLS.fold[k]   <- sum(test)*log(pred.MSE.OLS.fold) + length(m)*q*log(sum(test))
#     
#   }
#   
#   pred.MSE.GFI        <- mean((Y.pred.GFI - Y)^2)
#   pred.MAD.GFI        <- sum((rowMedians(abs(Y.pred.GFI-Y))^2))
#   pred.BIC.GFI        <- sum(pred.BIC.GFI.fold)
#   
#   pred.MSE.OLS        <- mean((Y.pred.OLS - Y)^2)
#   pred.MAD.OLS        <- sum((rowMedians(abs(Y.pred.OLS-Y))^2))
#   pred.BIC.OLS        <- sum(pred.BIC.OLS.fold)
#   
#   list("GFI"=c(pred.MSE.GFI, pred.MAD.GFI, pred.BIC.GFI), "OLS"=c(pred.MSE.OLS, pred.MAD.OLS, pred.BIC.OLS))
#   
# }





# admissible_subsets_CV <- function(Y, X, N = 100, steps = 2e2, burnin = 1e2, epsilonGrid, nFold = 5, 
#                                   Efficient=FALSE, useMCsample4h=FALSE, addoffset=FALSE, parallel=TRUE, n_cores=4){
#   
#   # print(getwd())
#   # source("SourceFunctions/Admissible_Subsets.R")
#   # source("SourceFunctions/Admissible_Subsets_Efficient.R")
#   require(foreach)
#   n <- ncol(Y); p <- nrow(X); q <- nrow(Y); 
#   fold            <- sample(1:nFold, n, replace = TRUE)
#   fold_prop       <- as.vector(table(fold))/n
#   
#   pred.MSE.OLS    <- matrix(NA, length(epsilonGrid), nFold)
#   pred.MSE.GFI    <- matrix(NA, length(epsilonGrid), nFold)
#   
#   pred.BIC.OLS    <- matrix(NA, length(epsilonGrid), nFold)
#   pred.BIC.GFI    <- matrix(NA, length(epsilonGrid), nFold)
#   
#   for (k in 1:nFold){
#     train         <- fold != k
#     test          <- fold == k
#     
#     Y.train       <- Y[,train]
#     X.train       <- X[,train]
#     
#     print(paste0("Fold: ", k))
#     
#     if (parallel){
#       
#       library(doParallel)
#       library(doSNOW)
#       cl                  <- makeSOCKcluster(n_cores)
#       registerDoSNOW(cl)
#       progress            <- function(nfin, tag) { cat(sprintf('tasks completed: %d; tag: %d\n', nfin, tag)) }
#       opts                <- list(progress=progress)
#       packages_req        <- c("matrixsampling", "tidyverse", "Rfast", "MSGLasso", "MASS")
#       
#       if (Efficient){
#         result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
#                                      .export = ls(globalenv()) ) %dopar% { admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                                            PropWeights=NULL, epsilon=epsi, true.Sigma = NULL, 
#                                                                            addoffset=addoffset, useMCsample4h=useMCsample4h)  }
#       } else{
#         result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
#                                      .export = ls(globalenv()) ) %dopar% { admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                                            PropWeights=NULL, epsilon=epsi, true.Sigma = NULL, 
#                                                                            addoffset=addoffset, useMCsample4h=useMCsample4h)  }
#       }
#       
#       parallel::stopCluster(cl)
#       
#     } else{
#       
#       if (Efficient){
#         result            <- foreach(epsi=epsilonGrid) %do% { admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                                     PropWeights=NULL, epsilon=epsi, true.Sigma = NULL, 
#                                                                     addoffset=addoffset, useMCsample4h=useMCsample4h)  }
#       } else{
#         result            <- foreach(epsi=epsilonGrid) %do% { admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                                    PropWeights=NULL, epsilon=epsi, true.Sigma = NULL, 
#                                                                    addoffset=addoffset, useMCsample4h=useMCsample4h)  }
#       }
#     }
#     
#     for (i in seq_along(epsilonGrid)){
#       
#       m                      <- result[[i]]$FinalModel
#       
#       print(paste0("For epsilon = ", epsilonGrid[i]))
#       cat("Chosen model:", m)
#       cat("\n")
#       
#       if (is.null(m)){
#         pred.MSE.GFI[i,k]    <-  pred.BIC.GFI[i,k]   <- NA ;
#         pred.MSE.OLS[i,k]    <-  pred.BIC.OLS[i,k]   <- NA ;
#       } else{
#         pred.GFI             <- Predict_MLM(Y.train=Y[, train], X.train=X[m, train, drop=FALSE], X.test=X[m, test, drop=FALSE], method="GFI", nsamp=1e3)$yhat
#         pred.OLS             <- Predict_MLM(Y.train=Y[, train], X.train=X[m, train, drop=FALSE], X.test=X[m, test, drop=FALSE], method="OLS", nsamp=1e3)$yhat
#         
#         pred.MSE.GFI[i,k]    <- mean((Y[, test] - pred.GFI)^2)
#         pred.MSE.OLS[i,k]    <- mean((Y[, test] - pred.OLS)^2)
#         
#         pred.BIC.GFI[i,k]    <- sum(test)*log(pred.MSE.GFI[i,k]) + length(m)*q*log(sum(test))
#         pred.BIC.OLS[i,k]    <- sum(test)*log(pred.MSE.OLS[i,k]) + length(m)*q*log(sum(test))
#         
#       }
#       
#     }
#     
#   }
#   
#   MSE.GFI             <- rowSums(pred.MSE.GFI %*% diag(fold_prop)) 
#   BIC.GFI             <- rowSums(pred.BIC.GFI %*% diag(fold_prop)) 
#   
#   MSE.OLS             <- rowSums(pred.MSE.OLS %*% diag(fold_prop)) 
#   BIC.OLS             <- rowSums(pred.BIC.OLS %*% diag(fold_prop)) 
#   
#   MSE.GFI.eps         <- epsilonGrid[which.min(MSE.GFI)]
#   BIC.GFI.eps         <- epsilonGrid[which.min(BIC.GFI)]
#   MSE.OLS.eps         <- epsilonGrid[which.min(MSE.OLS)]
#   BIC.OLS.eps         <- epsilonGrid[which.min(BIC.OLS)]
#   
#   optimal.eps         <- c(MSE.GFI.eps, BIC.GFI.eps, MSE.OLS.eps, BIC.OLS.eps)
#   names(optimal.eps)  <- c("MSE.GFI", "BIC.GFI", "MSE.OLS", "BIC.OLS")
#   
#   
#   list("optimal.eps" =  optimal.eps, "GFI"=cbind(MSE.GFI, BIC.GFI), "OLS"=cbind(MSE.OLS, BIC.OLS))
#   
# }







# admissible_subsets_CV <- function(Y, X, N = 100, steps = 2e2, burnin = 1e2, epsilonGrid, nFold = 5, 
#                                   Efficient=FALSE, useMCsample4h=FALSE, addoffset=FALSE, parallel=TRUE, n_cores=4, trueModel){
#   
#   
#   require(foreach); require(MASS);
#   n <- ncol(Y); p <- nrow(X); q <- nrow(Y); 
#   fold            <- sample(1:nFold, n, replace = TRUE)
#   fold_prop       <- as.vector(table(fold))/n
#   
#   pred.MSE.OLS    <- matrix(NA, length(epsilonGrid), nFold)
#   # pred.MSE.GFI    <- matrix(NA, length(epsilonGrid), nFold)
#   
#   pred.BIC.OLS    <- matrix(NA, length(epsilonGrid), nFold)
#   # pred.BIC.GFI    <- matrix(NA, length(epsilonGrid), nFold)
#   
#   MGLasso_Soln    <- quiet(InitialValue_MGLasso(Y.m=t(Y), X.m=t(X), lamG=seq(0.05,4,0.1), fold=5))
#   PropWeights     <- rowSums(MGLasso_Soln^2)
#   nnzweights      <- which(PropWeights != 0)
#   
#   for (k in 1:nFold){
#     train         <- fold != k
#     test          <- fold == k
#     
#     Y.train       <- Y[,train]
#     X.train       <- X[,train]
#     
#     print(paste0("Fold: ", k))
#     
#     if (parallel){
#       
#       library(doParallel)
#       library(doSNOW)
#       cl                  <- makeSOCKcluster(n_cores)
#       registerDoSNOW(cl)
#       progress            <- function(nfin, tag) { cat(sprintf('tasks completed: %d; tag: %d\n', nfin, tag)) }
#       opts                <- list(progress=progress)
#       packages_req        <- c("matrixsampling", "tidyverse", "Rfast", "MSGLasso", "MASS")
#       
#       if (Efficient){
#         result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
#                                      .export = ls(globalenv()) ) %dopar% { admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                                                                         PropWeights=PropWeights, epsilon=epsi, true.Sigma = NULL, 
#                                                                                                         addoffset=addoffset, useMCsample4h=useMCsample4h,
#                                                                                                         trueModel = trueModel)  }
#       } else{
#         result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
#                                      .export = ls(globalenv()) ) %dopar% { admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                                                               PropWeights=PropWeights, epsilon=epsi, true.Sigma = NULL, 
#                                                                                               addoffset=addoffset, useMCsample4h=useMCsample4h,
#                                                                                               trueModel = trueModel)  }
#       }
#       
#       parallel::stopCluster(cl)
#       
#     } else{
#       
#       if (Efficient){
#         result <- list()
#         for (i in seq_along(epsilonGrid)){
#           result[[i]] <- admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                       PropWeights=PropWeights, epsilon=epsilonGrid[i], true.Sigma = NULL, 
#                                                       addoffset=addoffset, useMCsample4h=useMCsample4h, trueModel = trueModel)
#         }
#         
#       } else{
#         result <- list()
#         for (i in seq_along(epsilonGrid)){
#           result[[i]] <- admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                       PropWeights=PropWeights, epsilon=epsilonGrid[i], true.Sigma = NULL, 
#                                                       addoffset=addoffset, useMCsample4h=useMCsample4h, trueModel = trueModel)
#         }
#         
#       }
#     }
#     
#     print(sum(test))
#     
#     for (i in seq_along(epsilonGrid)){
#       
#       m                      <- result[[i]]$FinalModel
#       all.m                  <- 
#       
#       print(paste0("For epsilon = ", epsilonGrid[i]))
#       cat("Chosen model:", m)
#       cat("\n")
#       
#       if (is.null(m)){
#         # pred.MSE.GFI[i,k]    <-  pred.BIC.GFI[i,k]   <- NA ;
#         pred.MSE.OLS[i,k]    <-  pred.BIC.OLS[i,k]   <- NA ;
#       } else{
#         
#         pred.OLS              <- Predict_MLM(Y.train=Y[, train, drop=FALSE], X.train=X[m, train, drop=FALSE],
#                                             X.test=X[m, test, drop=FALSE], method="OLS", nsamp=1e3)$yhat
#         # print("OLS is done")
#         # pred.GFI             <- Predict_MLM(Y.train=Y[, train, drop=FALSE], X.train=X[m, train, drop=FALSE], 
#                                                # X.test=X[m, test, drop=FALSE], method="GFI", nsamp=1e3)$yhat
#         
#         # pred.MSE.GFI[i,k]    <- mean((Y[, test, drop=FALSE] - pred.GFI)^2)
#         pred.MSE.OLS[i,k]    <- mean((Y[, test, drop=FALSE] - pred.OLS)^2)
#         
#         # pred.BIC.GFI[i,k]    <- sum(test)*log(pred.MSE.GFI[i,k]) + length(m)*log(sum(test))
#         pred.BIC.OLS[i,k]    <- sum(test)*log(pred.MSE.OLS[i,k]) + length(m)*log(sum(test))
#         
#       }
#       
#     }
#     
#   }
#   
#   MSE.GFI             <- rowSums(pred.MSE.GFI %*% diag(fold_prop)) 
#   BIC.GFI             <- rowSums(pred.BIC.GFI %*% diag(fold_prop)) 
#   
#   MSE.OLS             <- rowSums(pred.MSE.OLS %*% diag(fold_prop)) 
#   BIC.OLS             <- rowSums(pred.BIC.OLS %*% diag(fold_prop)) 
#   
#   MSE.GFI.eps         <- epsilonGrid[which.min(MSE.GFI)]
#   BIC.GFI.eps         <- epsilonGrid[which.min(BIC.GFI)]
#   MSE.OLS.eps         <- epsilonGrid[which.min(MSE.OLS)]
#   BIC.OLS.eps         <- epsilonGrid[which.min(BIC.OLS)]
#   
#   optimal.eps         <- c(MSE.GFI.eps, BIC.GFI.eps, MSE.OLS.eps, BIC.OLS.eps)
#   names(optimal.eps)  <- c("MSE.GFI", "BIC.GFI", "MSE.OLS", "BIC.OLS")
#   
#   
#   print(pred.MSE.OLS); 
#   print(MSE.OLS)
#   
#   
#   list("optimal.eps" =  optimal.eps, "GFI"=cbind(MSE.GFI, BIC.GFI), "OLS"=cbind(MSE.OLS, BIC.OLS))
#   
# }




# 
# admissible_subsets_CV <- function(Y, X, N = 100, steps = 2e2, burnin = 1e2, epsilonGrid, nFold = 5, 
#                                   Efficient=FALSE, useMCsample4h=FALSE, addoffset=FALSE, parallel=TRUE, n_cores=4, trueModel){
#   
#   
#   require(foreach); require(MASS);
#   n <- ncol(Y); p <- nrow(X); q <- nrow(Y); 
#   fold            <- sample(1:nFold, n, replace = TRUE)
#   fold_freq       <- table(fold)
#   
#   folds.occ       <- as.numeric(names(fold_freq))
#   
#   fold_prop       <- as.vector(fold_freq)/n
#   
#   nFold           <- length(folds.occ)
#   
#   print(paste0("Number of folds: ", nFold))
#   
#   print("Fold Proportions")
#   print(fold_prop)
#   
#   pred.MSE.OLS    <- matrix(NA, length(epsilonGrid), nFold)
#   pred.BIC.OLS    <- matrix(NA, length(epsilonGrid), nFold)
#   
#   pred.MSE.BMA    <- matrix(NA, length(epsilonGrid), nFold)
#   pred.BIC.BMA    <- matrix(NA, length(epsilonGrid), nFold)
#   
#   MGLasso_Soln    <- quiet(InitialValue_MGLasso(Y.m=t(Y), X.m=t(X), lamG=seq(0.05,4,0.1), fold=10))
#   PropWeights     <- rowSums(MGLasso_Soln^2)
#   nnzweights      <- which(PropWeights != 0)
#   
#   for (k in 1:nFold){
#     train         <- fold != folds.occ[k]
#     test          <- fold == folds.occ[k]
#     
#     Y.train       <- Y[,train]
#     X.train       <- X[,train]
#     
#     print(paste0("Fold: ", k))
#     print(paste0("FoldValue: ", folds.occ[k]))
#     
#     if (parallel){
#       
#       library(doParallel)
#       library(doSNOW)
#       cl                  <- makeSOCKcluster(n_cores)
#       registerDoSNOW(cl)
#       progress            <- function(nfin, tag) { cat(sprintf('tasks completed: %d; tag: %d\n', nfin, tag)) }
#       opts                <- list(progress=progress)
#       packages_req        <- c("matrixsampling", "tidyverse", "Rfast", "MSGLasso", "MASS")
#       
#       if (Efficient){
#         result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
#                                      .export = ls(globalenv()) ) %dopar% { admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                                                                         PropWeights=PropWeights, epsilon=epsi, true.Sigma = NULL, 
#                                                                                                         addoffset=addoffset, useMCsample4h=useMCsample4h,
#                                                                                                         trueModel = trueModel)  }
#       } else{
#         result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
#                                      .export = ls(globalenv()) ) %dopar% { admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                                                               PropWeights=PropWeights, epsilon=epsi, true.Sigma = NULL, 
#                                                                                               addoffset=addoffset, useMCsample4h=useMCsample4h,
#                                                                                               trueModel = trueModel)  }
#       }
#       
#       parallel::stopCluster(cl)
#       
#     } else{
#       
#       if (Efficient){
#         result <- list()
#         for (i in seq_along(epsilonGrid)){
#           result[[i]] <- admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                                       PropWeights=PropWeights, epsilon=epsilonGrid[i], true.Sigma = NULL, 
#                                                       addoffset=addoffset, useMCsample4h=useMCsample4h, trueModel = trueModel)
#         }
#         
#       } else{
#         result <- list()
#         for (i in seq_along(epsilonGrid)){
#           result[[i]] <- admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
#                                             PropWeights=PropWeights, epsilon=epsilonGrid[i], true.Sigma = NULL, 
#                                             addoffset=addoffset, useMCsample4h=useMCsample4h, trueModel = trueModel)
#         }
#         
#       }
#     }
#     
#     print(paste0("Number of test samples in fold: ", k))
#     print(sum(test))
#     
#     for (i in seq_along(epsilonGrid)){
#       
#       Opt.model              <- result[[i]]$Opt.model
#       Opt.B                  <- result[[i]]$Opt.B
#       Avg.B                  <- result[[i]]$Avg.B
#       Avg.size               <- result[[i]]$Avg.size
#         
#       print(paste0("For epsilon = ", epsilonGrid[i]))
#       cat("Chosen model:", Opt.model)
#       cat("\n")
#       
#       if (is.null(Opt.model)){
#         pred.MSE.BMA[i,k]    <-  pred.BIC.BMA[i,k]   <- NA ;
#         pred.MSE.OLS[i,k]    <-  pred.BIC.OLS[i,k]   <- NA ;
#       } else{
#         
#         pred.OLS             <- Opt.B %*% X[, test, drop=FALSE] 
#         pred.BMA             <- Avg.B %*% X[, test, drop=FALSE]
#         
#         pred.MSE.BMA[i,k]    <- mean((Y[, test, drop=FALSE] - pred.BMA)^2)
#         pred.MSE.OLS[i,k]    <- mean((Y[, test, drop=FALSE] - pred.OLS)^2)
#         
#         pred.BIC.BMA[i,k]    <- sum(test)*log(pred.MSE.BMA[i,k]) + Avg.size*log(sum(test))
#         pred.BIC.OLS[i,k]    <- sum(test)*log(pred.MSE.OLS[i,k]) + length(Opt.model)*log(sum(test))
#         
#       }
#       
#     }
#     
#   }
#   
#   
#   print("Fold calculations at a glance")
#   print(pred.MSE.BMA); print(pred.BIC.BMA);
#   print(pred.MSE.OLS); print(pred.BIC.OLS);
#   
#   MSE.BMA             <- rowSums(pred.MSE.BMA %*% diag(fold_prop)) 
#   BIC.BMA             <- rowSums(pred.BIC.BMA %*% diag(fold_prop)) 
#   
#   MSE.OLS             <- rowSums(pred.MSE.OLS %*% diag(fold_prop)) 
#   BIC.OLS             <- rowSums(pred.BIC.OLS %*% diag(fold_prop)) 
#   
#   MSE.BMA.eps         <- epsilonGrid[which.min(MSE.BMA)]
#   BIC.BMA.eps         <- epsilonGrid[which.min(BIC.BMA)]
#   MSE.OLS.eps         <- epsilonGrid[which.min(MSE.OLS)]
#   BIC.OLS.eps         <- epsilonGrid[which.min(BIC.OLS)]
#   
#   optimal.eps         <- c(MSE.BMA.eps, BIC.BMA.eps, MSE.OLS.eps, BIC.OLS.eps)
#   names(optimal.eps)  <- c("MSE.BMA", "BIC.BMA", "MSE.OLS", "BIC.OLS")
#   
#   
#   print(pred.MSE.OLS); 
#   print(MSE.OLS);
#   
#   
#   list("optimal.eps" =  optimal.eps, "BMA"=cbind(MSE.BMA, BIC.BMA), "OLS"=cbind(MSE.OLS, BIC.OLS))
#   
# }




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
  pred.BIC.OLS    <- matrix(NA, length(epsilonGrid), nFold)
  
  pred.MSE.BMA    <- matrix(NA, length(epsilonGrid), nFold)
  pred.BIC.BMA    <- matrix(NA, length(epsilonGrid), nFold)
  
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
                                     .export = ls(globalenv()) ) %dopar% { admissible_subsets_Efficient(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
                                                                                                        PropWeights=PropWeights, epsilon=epsi, true.Sigma = NULL, 
                                                                                                        addoffset=addoffset, useMCsample4h=useMCsample4h,
                                                                                                        trueModel = trueModel,MLasso=MLasso)  }
      } else{
        result            <- foreach(epsi=epsilonGrid, .options.snow=opts, .packages =packages_req, 
                                     .export = ls(globalenv()) ) %dopar% { admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
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
                                                      addoffset=addoffset, useMCsample4h=useMCsample4h, trueModel = trueModel, MLasso=MLasso)
        }
        
      } else{
        result <- list()
        for (i in seq_along(epsilonGrid)){
          result[[i]] <- admissible_subsets(Y=Y.train,X=X.train,N=N,steps=steps,burnin=burnin,
                                            PropWeights=PropWeights, epsilon=epsilonGrid[i], true.Sigma = NULL, 
                                            addoffset=addoffset, useMCsample4h=useMCsample4h, trueModel = trueModel, MLasso=MLasso)
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
      
      if (is.null(Opt.model)){
        pred.MSE.BMA[i,k]    <-  pred.BIC.BMA[i,k]   <- NA ;
        pred.MSE.OLS[i,k]    <-  pred.BIC.OLS[i,k]   <- NA ;
      } else{
        
        pred.OLS             <- Opt.B %*% X[, test, drop=FALSE] 
        pred.BMA             <- Avg.B %*% X[, test, drop=FALSE]
        
        pred.MSE.BMA[i,k]    <- mean((Y[, test, drop=FALSE] - pred.BMA)^2)
        pred.MSE.OLS[i,k]    <- mean((Y[, test, drop=FALSE] - pred.OLS)^2)
        
        # Error.comm.BMA       <- tcrossprod(Y[, test, drop=FALSE] - pred.BMA)
        # Error.comm.OLS       <- tcrossprod(Y[, test, drop=FALSE] - pred.OLS)
        
        pred.BIC.BMA[i,k]    <- sum(test)*log(pred.MSE.BMA[i,k]) + Avg.size*log(sum(test))
        pred.BIC.OLS[i,k]    <- sum(test)*log(pred.MSE.OLS[i,k]) + length(Opt.model)*log(sum(test))
        
        
        
      }
      
    }
    
  }
  
  
  print("Fold calculations at a glance")
  print(pred.MSE.BMA); print(pred.BIC.BMA);
  print(pred.MSE.OLS); print(pred.BIC.OLS);
  
  # MSE.BMA             <- rowSums(pred.MSE.BMA %*% diag(fold_prop)) 
  # BIC.BMA             <- rowSums(pred.BIC.BMA %*% diag(fold_prop)) 
  # 
  # MSE.OLS             <- rowSums(pred.MSE.OLS %*% diag(fold_prop)) 
  # BIC.OLS             <- rowSums(pred.BIC.OLS %*% diag(fold_prop)) 
  
  MSE.BMA             <- apply(pred.MSE.BMA, 1, weighted.mean, w=fold_prop, na.rm = FALSE)
  BIC.BMA             <- apply(pred.BIC.BMA, 1, weighted.mean, w=fold_prop, na.rm = FALSE)
  
  MSE.OLS             <- apply(pred.MSE.OLS, 1, weighted.mean, w=fold_prop, na.rm = FALSE)
  BIC.OLS             <- apply(pred.BIC.OLS, 1, weighted.mean, w=fold_prop, na.rm = FALSE)
  
  MSE.BMA.eps         <- epsilonGrid[which.min(MSE.BMA)]
  BIC.BMA.eps         <- epsilonGrid[which.min(BIC.BMA)]
  MSE.OLS.eps         <- epsilonGrid[which.min(MSE.OLS)]
  BIC.OLS.eps         <- epsilonGrid[which.min(BIC.OLS)]
  
  optimal.eps         <- c(MSE.BMA.eps, BIC.BMA.eps, MSE.OLS.eps, BIC.OLS.eps)
  names(optimal.eps)  <- c("MSE.BMA", "BIC.BMA", "MSE.OLS", "BIC.OLS")
  
  
  print(pred.MSE.OLS); 
  print(MSE.OLS);
  
  
  list("optimal.eps" =  optimal.eps, "BMA"=cbind(MSE.BMA, BIC.BMA), "OLS"=cbind(MSE.OLS, BIC.OLS))
  
}





