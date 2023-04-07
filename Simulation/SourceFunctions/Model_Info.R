

# Get the average value of the h-function by sampling from Matrix-T 
# distribution centered at the least squared estimates. Also compute the
# the generalized fiducial distribution of any particular model

model_info <- function(Y, X, L.X, active.set, epsilon, N=100, true.Sigma=NULL, useMCsample=FALSE){
  
  require(matrixsampling)
  n <- ncol(Y); p <- nrow(X); q <- nrow(Y); pm <- length(active.set);
  Xm            <- X[active.set, ,drop=FALSE]
  XmXmt         <- tcrossprod(Xm)
  XmYt          <- tcrossprod(Xm,Y)
  chol.Xm       <- chol(XmXmt)
  Bhat.m        <- t(backsolve(chol.Xm,forwardsolve(t(chol.Xm), XmYt)))
  # print(Bhat.m)
  
  if (is.null(true.Sigma)){
    Yhat          <- Y - Bhat.m%*%Xm
    Sighat.m      <- tcrossprod(Yhat)
  }else{
    Sighat.m      <- true.Sigma*n
  }
  
  
  # chol.Sig          <- chol(Sighat.m)
  
  svd.Sig    <- svd(Sighat.m)
  
  # print("Minimum eigen value of Sigmahat")
  # print(min(svd.Sig$d))
  # print("Condition number of Sigmahat")
  # print(max(svd.Sig$d)/min(svd.Sig$d))
  # print(min(diag(chol.Sig)))
  

  if (min(svd.Sig$d) < 1e-1)  {
    avg_h           <- 0
  } else{
    chol.Sig        <- chol(Sighat.m)
    Sighat.m.inv    <- chol2inv(chol.Sig)
    if (useMCsample) {
      nu            <- n - (pm + q) + 1
      sampled_B     <- rmatrixt(N, nu, Bhat.m, Sighat.m, solve(XmXmt), checkSymmetry=F)
      
      L             <- L.X/min(svd.Sig$d)
      # hcalc         <- lapply(1:N, function(index) {compute_h_function_MatrixForm(Bm=sampled_B[,,index], X=X, Sigma.Inv=Sighat.m.inv, L=L,
      #                                                                             active.set=active.set, eps=epsilon,
      #                                                                             maxit = 1e2, threshold=1e-5) } )

      hcalc         <- lapply(1:N, function(index) {compute_h_function_MatrixForm_Stable(Bm=sampled_B[,,index], X=X, Sig.half=chol.Sig, L=L,
                                                                                         active.set=active.set, eps=epsilon,
                                                                                         maxit = 5e2, threshold=1e-5) } )
      
      # print(mean(sapply(hcalc, function(x) x$obj)))
      # print(sd(sapply(hcalc, function(x) x$obj)))
      avg_h         <- mean(sapply(hcalc, function(x) x$hval))
      
    } else{
      L             <- L.X/min(svd.Sig$d)
      # hcalc         <- compute_h_function_MatrixForm(Bm=Bhat.m, X=X, Sigma.Inv=Sighat.m.inv, L=L, active.set=active.set,
      #                                                eps=epsilon, maxit = 1e2, threshold=1e-5)
      
      hcalc         <- compute_h_function_MatrixForm_Stable(Bm=Bhat.m, X=X, Sig.half=chol.Sig, L=L, active.set=active.set,
                                                            eps=epsilon, maxit = 5e2, threshold=1e-5)
      
      avg_h         <- hcalc$hval
      # print(hcalc$obj)
      # print(determinant(Sighat.m)$modulus)
    }
  }

  
  if (avg_h == 0){
    logf        <- -Inf
  } else{
    logAvg_h    <- log(avg_h)
    logf        <- pm*q*log(pi) + logAvg_h - 0.5*(n-pm-q)*(determinant(Sighat.m)$modulus) + 
                   0.25*q*(q-1)*log(pi) + sum(lgamma((n-pm-(1:q)+1)/2))
  }

  # c(avg_h, logf, determinant(Sighat.m)$modulus, pm*q*log(pi), sum(lgamma((n-pm-(1:q)+1)/2)),  0.5*(n-pm-q)*(determinant(Sighat.m)$modulus))
  c(avg_h, logf)
}
