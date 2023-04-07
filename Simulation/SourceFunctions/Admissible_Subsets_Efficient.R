
# The main function to walk through all possible choice of model by
# MCMC algoeithm, specified by the proposal models. This function 
# returns a MCMC chain with the selected model at each stage

admissible_subsets_Efficient <- function(Y,X,N=100,steps=15e3,burnin=5e3,PropWeights=NULL,epsilon, 
                                         true.Sigma=NULL, addoffset=TRUE, useMCsample4h=FALSE, MLasso=TRUE,
                                         trueModel){
  
  require(tidyverse)
  n <- ncol(Y); p <- nrow(X); q <- nrow(Y); 
  
  if (is.null(PropWeights)){
    if (MLasso){
      MLasso_Soln    <- InitialValue_MLasso(Y.m=t(Y), X.m=t(X))
      PropWeights    <- unname(rowSums(MLasso_Soln^2))
    } else{
      MGLasso_Soln   <- quiet(InitialValue_MGLasso(Y.m=t(Y), X.m=t(X), lamG=seq(0.05,10,0.1), fold=10))
      PropWeights    <- unname(rowSums(MGLasso_Soln^2))
    }
    if (addoffset){
      PropWeights  <- PropWeights + {n*p}^{-2}
    } 
    
  }
  
  nnzweights     <- which(PropWeights != 0)
  print(paste0("The length of initial model chosen: ", length(nnzweights)))
  
  ord.weights   <- order(PropWeights, decreasing=T)
  
  # print("Order of the weights")
  # print(ord.weights)
  
  # tweLarge      <- Rfast::nth(PropWeights,20,descending=T)
  # first20th     <- which(PropWeights >= tweLarge)
  
  # print("Printing 20th Large coefficients")
  # print(first20th)
  # print(trueModel)
  # Sys.sleep(4)
  
  XXt           <- tcrossprod(X)
  XYt           <- tcrossprod(X,Y)
  L.X           <- max(svd(XXt)$d)
  
  m             <- which.max(PropWeights)
  check.step    <- 1
  
  repeat{
    Info        <- model_info(Y = Y, X = X, L.X = L.X, active.set = m, epsilon=epsilon, N=N, 
                                true.Sigma=true.Sigma, useMCsample=useMCsample4h)
    if (Info[1]==1){
      print("Found an admissible model")
      print(m)
      # Sys.sleep(7)
      break
    }
    
    check.step   <- check.step + 1
    
    if (check.step >= n){
      print("Couldn't find any admissible model of size 1")
      
      Bhat.avg      <- matrix(0,q,p)
      Bhat.opt      <- matrix(0,q,p)
      Sighat.avg    <- (1/n)*(tcrossprod(Y - Bhat.avg %*% X))
      Sighat.opt    <- (1/n)*(tcrossprod(Y - Bhat.opt %*% X))
      
      det.Sighat.avg<- determinant(Sighat.avg)
      det.Sighat.opt<- determinant(Sighat.opt)
      
      BIC.avg       <- n*as.vector(det.Sighat.avg$modulus) + n*q*(log(2*pi)+1) + {q*(q+1)/2}*log(n)
      BIC.opt       <- n*as.vector(det.Sighat.opt$modulus) + n*q*(log(2*pi)+1) + {q*(q+1)/2}*log(n)
      
      return(list("FinalModel"=NULL, "AcceptRatio"=NULL, "ModelProb"=NULL, "lookup.dat"=NULL,
                  "Opt.B"=matrix(0,q,p), "Avg.B"=matrix(0,q,p), "Opt.model"= integer(0), "Avg.size" = 0,
                  "Avg.BIC"= BIC.avg, "Opt.BIC"= BIC.opt))
      break
    }
    
    # mnext        <- which(ord.weights == check.step)
    # mnext        <- which(PropWeights == tweLarge[check.step])
    mnext          <- ord.weights[check.step]
    # m            <- mnext
    m              <- sort(c(m, mnext))
    print("models updated")
    print(m)
    # Sys.sleep(3)
    
  }
  
  # m             <- which.max(colSums(true.B^2))
  
  # m             <- trueModel
  
  if ((length(setdiff(m, nnzweights)) > 0) & (addoffset==FALSE)){
    print("Initial chosen model is not a subset of non-zero weights")
    print("Adding offset to solve the issue of NA in probability")
    PropWeights  <- PropWeights + {n*p}^{-2}
    nnzweights     <- which(PropWeights != 0)
  }
  
  
  char.m        <- paste0(m, collapse = "|")
  cat("Initial model:", m)
  cat("\n")
  
  Info          <- model_info(Y = Y, X = X, L.X = L.X, active.set = m, epsilon=epsilon, N=N, 
                              true.Sigma=true.Sigma, useMCsample=useMCsample4h)
  oldlogf       <- Info[2]
  
  lookup.dat    <- data.frame(model=char.m, logf=Info[2], hfun=Info[1])
  
  chain         <- matrix(0, steps-burnin, p)
  AcceptRatio   <- 0
  
  for(t in 1:steps){
    
    if (t %% 100 == 0){
      print(paste(t, "iteration done"))
    }
    
    # print(paste("Updates for step", t))
    
    m_old       <- m
    
    # cat("Old model:", m_old, sep="\n")
    
    # Propose a new model
    AddVar      <- sample(c("Add", "Remove", "Exchange"),1)
    
    # cat("Requested to", AddVar, sep="\n")
    
    if (length(setdiff(nnzweights, m_old))==0 & (AddVar != "Remove")){
      m           <- m_old
    } else{
      
      if (any(PropWeights[m_old] == 0)){
        print(m_old)
        print(setdiff(m_old, nnzweights))
      }
      
      NewModel    <- proposal(addvar = AddVar, cur.model = m_old, n_covar = p, weights = PropWeights)
      m           <- NewModel$newmodel
      # print("Is the selected model subset of true model")
      # print(all(m %in% trueModel))
      # print("printing the size of the current model")
      # print(length(m))
      
    }
    
    # cat("Proposed model", m, sep="\n")
    
    char.m      <- paste0(m, collapse = "|")
    
    looked.logf <- lookup.dat %>% dplyr::filter(model==char.m) %>% pull(logf)
    looked.hfun <- lookup.dat %>% dplyr::filter(model==char.m) %>% pull(hfun)
    
    # cat("After looking up from the table, logf is", looked.logf, sep="\n")
    
    if (isTRUE(all.equal(m, m_old))){
        # Same model
      # print("Same model")
      newlogf    <- oldlogf
    } else{
      if (is_empty(looked.logf)){
        # Need computation
        # print("Not visited before : needs computation")
        Info       <- model_info(Y = Y, X = X, L.X = L.X, active.set = m, epsilon=epsilon, N=N, 
                                 true.Sigma=true.Sigma, useMCsample=useMCsample4h)
        newlogf    <- Info[2]
        newhfun    <- Info[1]
        # print(newhfun)
        # Updated the look up table
        lookup.dat <- lookup.dat %>% add_row(model=char.m, logf=Info[2], hfun=Info[1])
      }else{
        # Already visited
        # print("Visited before : No need of computation")
        newlogf    <- looked.logf
        newhfun    <- looked.hfun
      }
    }
    
    
    # cat("h function for the proposed model", newhfun, sep="\n")
    
    logrho      <-  newlogf - oldlogf + log(NewModel$qratio)
    
    # print(logrho)
    
    # cat("newlogf :", newlogf, sep="\n")
    # cat("oldlogf :", oldlogf, sep="\n")
    # cat("logrho :", logrho, sep="\n")
    
    if (log(runif(1)) < logrho){
      oldlogf       <- newlogf
      # print("****************************************")
      # print("ACCEPTED")
      # print("****************************************")
      if (t > burnin){
        AcceptRatio <- AcceptRatio + 1
      }
    } else{
      m             <- m_old
      # print("****************************************")
      # print("REJECTED")
      # print("****************************************")
    }
    
    # Sys.sleep(7)
    
    # cat("Updated model:", m)
    # cat("\n")
    
    if (t > burnin){
      chain[t-burnin,m] <- 1
    }
  }
  
  models        <- lapply(1:(steps-burnin), function(x) which(chain[x,]==1))
  models.char   <- unlist(lapply(models, paste, collapse = " "))
  m.final       <- models[[which(models.char == names(which.max(table(models.char))))[1]]]
  
  mchar.freq    <- table(models.char)
  mchar.all     <- names(mchar.freq)
  mchar.len     <- stringr::str_count(mchar.all, "\\w+")
  mchar.opt.ind <- which.max(mchar.freq)
  
  m.all         <- lapply(seq_along(mchar.all), function(i) as.numeric(stringr::word(mchar.all[i], 1:mchar.len[i])) )
  m.opt         <- m.all[[mchar.opt.ind]]
  
  Bhat.all      <- lapply(m.all, function(model) {X.sub         <- X[model, , drop=FALSE];
  fit           <- lm(t(Y) ~ -1 + t(X.sub) ); 
  Bhat          <- matrix(0, nrow(Y), nrow(X)) ; 
  Bhat[, model] <- t(unname(fit$coef)) ;
  Bhat
  })
  
  mchar.prob    <- as.vector(mchar.freq/sum(mchar.freq))
  Bhat.avg      <- Reduce("+", lapply(seq_along(Bhat.all), function(i) Bhat.all[[i]]*mchar.prob[i]) )
  Bhat.opt      <- Bhat.all[[mchar.opt.ind]]
  
  avg.size      <- mean(mchar.len)
  
  Sighat.avg    <- (1/(n-avg.size))*(tcrossprod(Y - Bhat.avg %*% X))
  Sighat.opt    <- (1/(n-length(m.opt)))*(tcrossprod(Y - Bhat.opt %*% X))
  
  BIC.avg       <- n*as.vector(determinant(Sighat.avg)$modulus) + n*q*(log(2*pi)+1) + {q*mean(mchar.len) + q*(q+1)/2}*log(n)
  BIC.opt       <- n*as.vector(determinant(Sighat.opt)$modulus) + n*q*(log(2*pi)+1) + {q*length(m.opt) + q*(q+1)/2}*log(n)
  
  AcceptRatio   <- AcceptRatio/steps
  
  
  list("FinalModel"=m.final, "AcceptRatio"=AcceptRatio, "ModelFreq"=table(models.char), "lookup.dat"=lookup.dat,
       "Opt.B"=Bhat.opt, "Avg.B"=Bhat.avg, "Opt.model"= m.opt, "Avg.size" = avg.size, 
       "Avg.BIC"= BIC.avg, "Opt.BIC"= BIC.opt)
  
}
