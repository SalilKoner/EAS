
cross_validation <- function(Y, X, model, nfold){
  n <- ncol(Y); p <- nrow(X); q <- nrow(Y); 
  K               <- nfold
  #set.seed(919)
  fold            <- sample(1:K,n,replace=TRUE)
  table(fold)
  Y.hat           <- matrix(NA,q,n)
  
  for(k in 1:K){
    train         <- fold != k
    test          <- fold == k
    X.train       <- X[model, train, drop=FALSE]
    Y.train       <- Y[,      train, drop=FALSE]
    X.test        <- X[model, test,  drop=FALSE]
    Y.test        <- Y[,      test,  drop=FALSE]
    fit           <- lm(t(Y.train) ~ -1 + t(X.train))
    Y.hat[, test] <- t(fit$coefficients) %*% X.test     
  }
  
  MSE             <- mean((Y-Y.hat)^2)
  MAD             <- median(abs(Y-Y.hat))
  
  c("MSE" = MSE, "MAD" = MAD)
  
}
