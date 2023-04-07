

# Computing the value of the h-function to determine the 
# epsilon-admissibility of a particular model through Projected Gradient Descent

compute_h_function_PGD <- function(B, Q, R, S, L, eps, active.set, q, p,threshold=1e-4, maxit = 5e2){
  
  if (length(active.set)==1){
    obj.val     <- sqrt(0.5*sum(as.vector(B)*as.vector(S %*% as.vector(B))))
    nnz         <- NULL
    exit.code   <- 0
    flips       <- NULL
  }
  
  else{
    require(Rfast)
    c             <- as.vector(crossprod(R, as.vector(B)))
    nm            <- which.min(colSums(B^2))
    nnz           <- active.set[-nm]
    nnp           <- as.vector(outer(1:q, (nnz-1)*q, FUN="+"))
    b.sig         <- as.vector(B[,-nm]) 
    Qb            <- as.vector(.subset(Q, 1:(q*p), nnp) %*% b.sig)
    btQb          <- 0.5*sum(Qb[nnp]*b.sig) 
    ctb           <- sum(c[nnp]*b.sig) 
    obj.cur       <- btQb - ctb
    const         <- 0.5*sum(as.vector(B)*as.vector(S %*% as.vector(B)))
    K             <- length(active.set)-1
    flips         <- NULL
    exit.code     <- 0
    print("Initially Objective")
    print(obj.cur + const)
    
    for (iter in 1:maxit){
      obj.prev    <- obj.cur
      b           <- - (Qb - c)/L
      b[nnp]      <- b[nnp] + b.sig
      colnorm2    <- colSums(matrix(b^2,q,p))
      # print(colnorm2)
      kth.large   <- Rfast::nth(colnorm2,K,descending=T)
      nnz         <- which(colnorm2 >= kth.large)
      nnp         <- as.vector(outer(1:q, (nnz-1)*q, FUN="+"))
      b.sig       <- b[nnp]           
      Qb          <- as.vector(.subset(Q, 1:(q*p), nnp) %*% b.sig)
      btQb        <- 0.5*sum(Qb[nnp]*b.sig) 
      ctb         <- sum(c[nnp]*b.sig) 
      obj.cur     <- btQb - ctb
      obj.diff    <- obj.prev - obj.cur
      obj.val     <- sqrt(obj.cur + const)
      if (obj.diff < 0){
        flips     <- c(flips, iter)
      }
      else if (obj.diff < threshold){
        # print(paste0('converged at iteration ', iter))
        exit.code  <- 0
        break
      }  
      print(iter)
    }
  }
  list("hval"=as.integer(obj.val > eps), "nnz"=nnz, "flips"=flips, "ex.code"=exit.code, "obj"=obj.val)
}




PGD_naive <- function(B.M, Q, S, c, active.set, eps, maxit = 5e2, L = 1, threshold=1e-7){
  require(Rfast); require(Matrix);
  n             <- nrow(B.M) # dimension of the multivariate response Y
  p             <- as.integer(nrow(Q)/n) # number of covariates in the full model
  nm            <- which.min(colSums(B.M^2))
  nnz           <- active.set[-nm]
  #print(nnz)
  B.cur         <- matrix(0,n,p)
  B.cur[,nnz]   <- B.M[,-nm]
  Qb            <- as.vector(Q%*%as.vector(B.cur))
  #print(Qb)
  btQb          <- 0.5*sum(Qb*as.vector(B.cur))
  #print(paste("Qcomp=", btQb))
  ctb           <- sum(c*as.vector(B.cur))
  #print(paste("Lcomp=", ctb))
  obj.cur       <- btQb - ctb
  const         <- 0.5*sum(as.vector(B.M)*as.vector(S %*% as.vector(B.M)))
  K             <- length(active.set)-1
  flips         <- NULL
  #print(obj.cur)
  
  for (iter in 1:maxit){
    obj.prev    <- obj.cur
    B.cur       <- B.cur - matrix((Qb - c)/L, n, p)
    colnorm2    <- colSums(B.cur^2)
    kth.large   <- Rfast::nth(colnorm2,K,descending=T)
    #print(kth.large)
    nnz         <- which(colnorm2 >= kth.large)
    B.cur[,-nnz]<- 0
    Qb          <- as.vector(Q%*%as.vector(B.cur))
    btQb        <- 0.5*sum(Qb*as.vector(B.cur))
    ctb         <- sum(c*as.vector(B.cur))
    obj.cur     <- btQb - ctb
    obj.diff    <- obj.prev - obj.cur
    obj.val     <- sqrt(obj.cur + const)
    if (obj.diff < 0){
      flips     <- c(flips, iter)
    }
    else if (obj.diff < threshold){
      exit.code  <- 0
      break
    }  
  }
  list("hval"=as.integer(obj.val > eps), "nnz"=nnz, "flips"=flips, "ex.code"=exit.code, "obj"=obj.val)
}  



# Computing the h function using the gradient descent in terms of the matrix form

compute_h_function_MatrixForm <- function(Bm, X, Sigma.Inv, L, active.set, eps, maxit = 5e2, threshold=1e-7){
  
  p         <- nrow(X) ; n  <- ncol(X) ; K <- length(active.set)-1; q <-nrow(Sigma.Inv);
  Xm        <- X[active.set, ,drop=FALSE]
  BmXm      <- Bm %*% Xm
  BmXmX     <- (Sigma.Inv %*% BmXm) %*% t(X)
  
  if (length(active.set) == 1){
    obj         <- 0.5*sum(diag(tcrossprod(BmXm)  %*% Sigma.Inv))
    nnz         <- NULL
    exit.code   <- "trivial"
    # print("trivial case")
    return(list("hval"=as.integer(obj > eps), "nnz"=nnz, "ex.code"=exit.code, "obj"=obj))
  }
  
  min.norm            <- active.set[which.min(colSums(Bm^2))] ; 
  B.cur               <- matrix(0, q, p) ; 
  B.cur[, active.set] <- Bm ; B.cur[, min.norm] <- 0;
  nnz                 <- setdiff(active.set, min.norm)
  
  obj                 <- 0.5*sum(diag(tcrossprod(BmXm - B.cur[, nnz, drop=FALSE] %*% X[nnz, , drop=FALSE])  %*% Sigma.Inv))
  exit.code           <- "maxIter"
  
  for (iter in 1:maxit){
    # print(iter)
    objold       <- obj
    grad         <- tcrossprod(Sigma.Inv %*% (B.cur[, nnz, drop=FALSE] %*% X[nnz, , drop=FALSE]), X) - BmXmX
    B.cur        <- B.cur - grad/L
    colnorm2     <- colSums(B.cur^2)
    kth.large    <- Rfast::nth(colnorm2,K,descending=T)
    nnz          <- which(colnorm2 >= kth.large)
    B.cur[,-nnz] <- 0
    obj          <- 0.5*sum(diag(tcrossprod(BmXm - B.cur[, nnz, drop=FALSE] %*% X[nnz, , drop=FALSE]) %*% Sigma.Inv ))
    objdiff      <- abs(obj - objold)
    # print(objdiff/objold); print(objdiff); print(nnz);
    if (objdiff/objold < threshold){
      diff.BX    <- BmXm - B.cur[, nnz, drop=FALSE] %*% X[nnz, ,drop=FALSE]
      obj        <- 0.5*sum(diag(tcrossprod(diff.BX)  %*% Sigma.Inv))
      exit.code  <- "conv"
      # print("converged")
      break
    }
    # print(iter)
  }
  list("hval"=as.integer(obj > eps), "nnz"=nnz, "ex.code"=exit.code, "obj"=obj)
}



compute_h_function_MatrixForm_Stable <- function(Bm, X, Sig.half, L, active.set, eps, maxit = 5e2, threshold=1e-7){
  
  p         <- nrow(X) ; n  <- ncol(X) ; K <- length(active.set)-1; q <-nrow(Sig.half);
  Xm        <- X[active.set, ,drop=FALSE]
  BmXm      <- Bm %*% Xm
  BmXmX     <- backsolve(Sig.half, forwardsolve(t(Sig.half), BmXm)) %*% t(X)
  
  if (length(active.set) == 1){
    obj         <- 0.5*sum(diag( backsolve(Sig.half, forwardsolve(t(Sig.half), tcrossprod(BmXm)))   ))
    nnz         <- NULL
    exit.code   <- "trivial"
    # print("trivial case")
    return(list("hval"=as.integer(obj > eps), "nnz"=nnz, "ex.code"=exit.code, "obj"=obj))
  }
  
  min.norm            <- active.set[which.min(colSums(Bm^2))] ; 
  B.cur               <- matrix(0, q, p) ; 
  B.cur[, active.set] <- Bm ; B.cur[, min.norm] <- 0;
  nnz                 <- setdiff(active.set, min.norm)
  
  BX                  <- B.cur[, nnz, drop=FALSE] %*% X[nnz, , drop=FALSE]
  obj                 <- 0.5*sum(diag( backsolve(Sig.half, forwardsolve(t(Sig.half), tcrossprod(BmXm - BX)))   ))
  exit.code           <- "maxIter"
  
  for (iter in 1:maxit){
    # print(iter)
    objold       <- obj
    grad         <- tcrossprod(backsolve(Sig.half, forwardsolve(t(Sig.half), BX)), X) - BmXmX
    B.cur        <- B.cur - grad/L
    colnorm2     <- colSums(B.cur^2)
    kth.large    <- Rfast::nth(colnorm2,K,descending=T)
    nnz          <- which(colnorm2 >= kth.large)
    B.cur[,-nnz] <- 0
    BX           <- B.cur[, nnz, drop=FALSE] %*% X[nnz, , drop=FALSE]
    obj          <- 0.5*sum(diag( backsolve(Sig.half, forwardsolve(t(Sig.half), tcrossprod(BmXm - BX)))   ))
    objdiff      <- abs(obj - objold)
    # print(objdiff/objold); print(objdiff); print(nnz);
    if (objdiff/objold < threshold){
      diff.BX    <- BmXm - BX
      obj        <- 0.5*sum(diag( backsolve(Sig.half, forwardsolve(t(Sig.half), tcrossprod(diff.BX)))   ))
      exit.code  <- "conv"
      # print("converged")
      break
    }
    # print(iter)
  }
  list("hval"=as.integer(obj > eps), "nnz"=nnz, "ex.code"=exit.code, "obj"=obj)
}





