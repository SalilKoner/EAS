
#  Get the Multivariate Group Lasso Estimate as a warm start 
#  to the EAS methodology for Grouped Variable Selection 

InitialValue_MGLasso <- function(Y.m, X.m, lamG, fold=5){
  
  require(MSGLasso)
  P              <- ncol(X.m) 
  Q              <- ncol(Y.m) 
  G              <- P
  R              <- 1
  gmax           <- 1
  cmax           <- Q
  
  GarrStarts     <- 0:(P-1)
  GarrEnds       <- 0:(P-1)
  RarrStarts     <- 0
  RarrEnds       <- Q-1
  
  tmp            <- FindingPQGrps(P, Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  PQ.grps        <- tmp$PQgrps
  
  tmp1           <- Cal_grpWTs(P, Q, G, R, gmax, PQ.grps)
  grpWTs         <- tmp1$grpWTs
  tmp2           <- FindingGRGrps(P, Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  GRgrps         <- tmp2$GRgrps
  
  # Pen_L          <- matrix(rep(0,P*Q),P,Q, byrow=TRUE)
  # Pen_G          <- matrix(rep(1,G*R),G,R, byrow=TRUE)
  
  assign("Pen_L", matrix(rep(0,P*Q),P,Q, byrow=TRUE), envir=globalenv())
  assign("Pen_G", matrix(rep(1,G*R),G,R, byrow=TRUE), envir=globalenv())
  
  # print(Pen_L)
  
  lamG.v         <- lamG
  
  try.cv         <- MSGLasso.cv(X.m, Y.m, grpWTs, Pen_L, Pen_G, PQ.grps, GRgrps, 0, lamG.v, fold=5, seed=1)
  
  # print("done here")
  
  grp_Norm0      <- matrix(rep(0, G*R), nrow=G, byrow=TRUE)
  MSGLassolam1   <- try.cv$lams.c[which.min(as.vector(try.cv$rss.cv))][[1]]$lam1
  MSGLassolamG   <- try.cv$lams.c[which.min(as.vector(try.cv$rss.cv))][[1]]$lam3
  MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)
  
  MGLasso.cv     <- MSGLasso(X.m, Y.m, grpWTs, Pen_L, Pen_G, PQ.grps, GRgrps, grp_Norm0,
                             MSGLassolam1, MSGLassolamG.m, Beta0=NULL)
  
  print(MSGLassolamG)
  
  MGLasso.cv$Beta
  
}