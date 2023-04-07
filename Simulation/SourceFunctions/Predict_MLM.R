
# This function provides the least square estimator of B and the 
# covariance matrix Sigma for tall data. 

Predict_MLM <- function(Y.train, X.train, X.test, method="GFI", nsamp=NULL, burn=NULL){
  
  require(MASS); require(matrixsampling);
  n.train <- ncol(Y.train); n.test <- ncol(X.test)
  p <- nrow(X.train); q <- nrow(Y.train);
  
  XmXmt                <- tcrossprod(X.train)
  XmYt                 <- tcrossprod(X.train,Y.train)
  chol.Xm              <- chol(XmXmt)
  Bhat.m               <- t(backsolve(chol.Xm,forwardsolve(t(chol.Xm), XmYt)))
  # Sighat.m             <- tcrossprod(Y.train - Bhat.m%*%X.train)
  # V.m                  <- Sighat.m/(n.train-p)
  # Y.test.pred          <- Bhat.m%*%X.test + t(matrix(mvrnorm(n.test, rep(0,q), V.m ), nrow=n.test, ncol=q ))
  Y.test.pred          <- Bhat.m%*%X.test
  
  if (method=="GFI"){
    nsamp              <- ifelse(is.null(nsamp), 1e2, nsamp)
    burn               <- ifelse(is.null(burn),  5e3, burn)
    nu                 <- n.train - (p + q) + 1
    sample.B           <- rmatrixt(nsamp, nu, Bhat.m, Sighat.m, solve(XmXmt), checkSymmetry=F)
    sample.Sigma       <- sapply(1:nsamp, function(index) drop(rinvwishart(1, n.train, tcrossprod(Y.train - sample.B[,,index]%*%X.train ))),
                                 simplify = "array")

    samp.Y.test.pred   <- sapply(1:nsamp, function(index) sample.B[,,index]%*%X.test + 
                                                          t(matrix(mvrnorm(n.test, rep(0,q), sample.Sigma[,,index] ), nrow=n.test, ncol=q )), simplify = "array" )
    Y.test.pred        <- apply(samp.Y.test.pred, c(1,2), mean)
    Bhat.m             <- apply(sample.B, c(1,2), mean)
    V.m                <- apply(sample.Sigma, c(1,2), mean)
    
  }
  
  list("yhat"=Y.test.pred, "Bhat" = Bhat.m)
  
  # list("yhat"=Y.test.pred, "Bhat" = Bhat.m, "Vhat"=V.m)
  # list("Bhat" = Bhat.m, "Vhat"=V.m)
  
}

# Predict_MLM(Y.train=Y, X.train=Xm, X.test=Xm, method="OLS", nsamp=1e4)
     
