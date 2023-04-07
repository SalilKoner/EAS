

InitialValue_MLasso <- function(Y.m, X.m){
  require(glmnet)
  q        <- ncol(Y.m)
  cv.net   <- cv.glmnet(X.m, Y.m, family="mgaussian", type.measure = "mse", intercept = FALSE, standardize = FALSE)
  obj      <- glmnet(X.m, Y.m, family="mgaussian", type.measure = "mse", intercept = FALSE, lambda=cv.net$lambda.min, standardize = FALSE)
  beta     <- coef(obj)
  beta.mat <- as.matrix(do.call(cbind,  lapply(1:q, function(i) beta[[i]])))
  beta.mat[-1, ]
}

# cv.net   <- cv.glmnet(t(X), t(Y), family="mgaussian", type.measure = "mse", intercept = FALSE, standardize = FALSE)
# obj      <- glmnet(t(X), t(Y), family="mgaussian", intercept = FALSE, lambda=cv.net$lambda.min, standardize = FALSE)
# beta     <-  coef(obj)
# beta.mat <- as.matrix(do.call(cbind,  lapply(1:q, function(i) beta[[i]])))
# beta.mat <- as.matrix(cbind(beta[[1]], beta[[2]], beta[[3]], beta[[4]], beta[[5]])[-1,])
# which.max(rowSums(beta.mat^2))
# sum(rowSums(beta.mat^2) != 0)
