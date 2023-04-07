
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


get_summary <- function(Y.in, X.in, Y.out, X.out, B.hat, B.true, model.hat, model.true, estModelFreq = NULL){
  n <- ncol(Y.in) ; q <- nrow(Y.in) ; p <- nrow(X.in)
  
  true.model.char      <- paste(model.true, collapse = " ")
  
  if(!is.null(estModelFreq)){
    mchar.all          <- names(estModelFreq)
    prob.true.model    <- unname(estModelFreq[mchar.all ==  true.model.char]/sum(estModelFreq))
    if (length(prob.true.model) == 0){
      prob.true.model  <- 0
    }
  } else{
    prob.true.model    <- NA
  }
  
  is.MAPmodel.true     <- as.integer(isTRUE(all.equal(model.hat, model.true)))
  
  mod                  <- rep(0,p)
  mod[model.true]      <- 1
  
  mod.est              <- integer(p)
  mod.est[model.hat]   <- 1
  
  TP                   <- sum(mod == 1 & mod.est==1)
  FP                   <- sum(mod == 0 & mod.est==1)
  FN                   <- sum(mod == 1 & mod.est==0) 
  TN                   <- sum(mod == 0 & mod.est==0) 
  
  FDR <- FP/(TP+FP); FNR <- FN/(FN + TN); MP <- (FP+FN)/(p*q);
  
  SE.B                 <- sum(diag(tcrossprod(B.true - B.hat)))/(p*q)
  Y.in.hat            <- B.hat %*% X.in
  SE.in.est           <- sum(diag(tcrossprod(B.true %*% X.in - Y.in.hat)))/(n*q)
  SE.in.pred          <- sum(diag(tcrossprod(Y.in - Y.in.hat)))/(n*q)
  
  
  Y.out.hat           <- B.hat %*% X.out
  SE.out.est          <- sum(diag(tcrossprod(B.true %*% X.out - Y.out.hat)))/(n*q)
  SE.out.pred         <- sum(diag(tcrossprod(Y.out - Y.out.hat)))/(n*q)
  
  measure             <- c(SE.B, SE.in.est, SE.in.pred, SE.out.est, SE.out.pred, FDR,FNR,MP, prob.true.model, is.MAPmodel.true)
  
  unname(measure)
  
}


get_summary_MSRL <- function(Y.out, Y.out.hat, B.hat, model.hat, model.true, estModelFreq = NULL){
  
  n <- ncol(Y.out) ; q <- nrow(Y.out) ; p <- ncol(B.hat)
  
  true.model.char      <- paste(model.true, collapse = " ")
  
  if(!is.null(estModelFreq)){
    mchar.all          <- names(estModelFreq)
    prob.true.model    <- unname(estModelFreq[mchar.all ==  true.model.char]/sum(estModelFreq))
    if (length(prob.true.model) == 0){
      prob.true.model  <- 0
    }
  } else{
    prob.true.model    <- NA
  }
  
  is.MAPmodel.true     <- as.integer(isTRUE(all.equal(model.hat, model.true)))
  
  mod                  <- rep(0,p)
  mod[model.true]      <- 1
  
  mod.est              <- integer(p)
  mod.est[model.hat]   <- 1
  
  TP                   <- sum(mod == 1 & mod.est==1)
  FP                   <- sum(mod == 0 & mod.est==1)
  FN                   <- sum(mod == 1 & mod.est==0) 
  TN                   <- sum(mod == 0 & mod.est==0) 
  
  FDR                 <- FP/(TP+FP); FNR <- FN/(FN + TN); MP <- (FP+FN)/(p*q);
  
  SE.out.pred         <- sum(diag(tcrossprod(Y.out - Y.out.hat)))/(n*q)
  
  measure             <- c(NA, NA, NA, NA, SE.out.pred, FDR,FNR,MP, prob.true.model, is.MAPmodel.true)
  
  unname(measure)
  
}