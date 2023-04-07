args     <- commandArgs(TRUE)
dsn      <- as.integer(args[1])
n_cores  <- as.integer(args[2])
Exper    <- as.integer(args[3])
eps.st   <- as.numeric(args[4])
eps.ed   <- as.numeric(args[5])
eps.ln   <- as.integer(args[6])
fold.num <- as.integer(args[7])

print(Exper)

# dsn <- 16; Exper <- 4; n_cores <- 8; eps.st <- 0.1 ; eps.ed <- 5 ; eps.ln <- 24; fold.num <- 10;
# dsn <- 49; Exper <- 6; n_cores <- 8; eps.st <- 0.2 ; eps.ed <- 0.5 ; eps.ln <- 20; fold.num <- 10;


# Fixing the seed to ensure that we get the same data at each node
set.seed(seed=(09112001+dsn))

source("setup_BaiGhosh2018.R")

# Generate data
dat <- generate_data(isim=dsn, n=n, q=q, p=p, m0=m0, Sig.X=Gamma.mat, Sig.E=Sigma.mat)


source("SourceFunctions/Proposal_Model.R")
source("SourceFunctions/Model_Info.R")
source("SourceFunctions/Compute_h.R")
source("SourceFunctions/MGLasso_Estimate.R")
source("SourceFunctions/MLasso_Estimate.R")
source("SourceFunctions/Admissible_Subsets.R")
source("SourceFunctions/Admissible_Subsets_Efficient.R")
source("SourceFunctions/Admissible_Subsets_CV.R")
source("SourceFunctions/Predict_MLM.R")
source("SourceFunctions/utilities.R")


Y                   <- dat$Y ; 
X                   <- dat$X ; 
true.model          <- dat$true.active ; 
true.B              <- dat$true.B; 
eps.grid            <- seq(eps.st, eps.ed, length.out = eps.ln)

print(eps.grid) ; print(fold.num);

st.CVEAS            <- proc.time()
# EAS method with epsilon chosen through 10-fold cross validation
cvfold              <- admissible_subsets_CV(Y=Y, X=X, N = 100, steps = 500, burnin = 200, epsilonGrid=eps.grid, 
                                             nFold = fold.num, Efficient=TRUE, useMCsample4h=FALSE, addoffset=FALSE, 
                                             parallel=TRUE, n_cores = n_cores, trueModel=true.model, MLasso=TRUE)

# Optimal chosen epsilon
print(cvfold$optimal.eps)

cl                  <- makeSOCKcluster(length(cvfold$optimal.eps))
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed: %d; tag: %d\n', nfin, tag)) }
opts                <- list(progress=progress)
packages_req        <- c("matrixsampling", "tidyverse", "Rfast", "MSGLasso", "MASS", "glmnet")
final.model.CV      <- foreach(epsi=cvfold$optimal.eps, .options.snow=opts, .packages =packages_req) %dopar% {
                                         admissible_subsets_Efficient(Y=Y, X=X, N = 1000, steps = 1e4, burnin = 5e3, PropWeights = NULL,
                                                                      epsilon = epsi, true.Sigma = NULL, addoffset=FALSE, 
                                                                      useMCsample4h = FALSE, MLasso=TRUE) } 

parallel::stopCluster(cl)
time.CVEAS          <- proc.time()[3]-st.CVEAS[3]
print("Done EAS-CV")

# EAS method with epsilon chosen through BIC
library(foreach); library(doSNOW) ; library(doParallel);

st.BICEAS           <- proc.time()
cl                  <- makeSOCKcluster(n_cores)
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed: %d; tag: %d\n', nfin, tag)) }
opts                <- list(progress=progress)
packages_req        <- c("matrixsampling", "tidyverse", "Rfast", "MSGLasso", "glmnet", "MASS")
BIC.models          <- foreach(epsi=eps.grid, .options.snow=opts, .packages =packages_req) %dopar% {
                                      admissible_subsets_Efficient(Y=Y,X=X,N=100,steps=5e3,burnin=2e3, PropWeights=NULL,
                                                                   epsilon=epsi, true.Sigma=NULL, addoffset=FALSE,
                                                                   useMCsample4h=FALSE, MLasso=TRUE) }

stopCluster(cl)

BIC.est              <- t(sapply(seq_along(BIC.models), function(i) c(BIC.models[[i]]$Opt.BIC, 
                                                                      BIC.models[[i]]$Avg.BIC, 
                                                                      length(BIC.models[[i]]$Opt.model) )) )
print(BIC.est)
min.BIC              <- min(BIC.est[,1])
min.BIC.index        <- which(BIC.est[,1] == min.BIC)
eps.chosen           <- eps.grid[min.BIC.index]
print(eps.chosen)

final.model.BIC      <- BIC.models[[sample(min.BIC.index, 1)]]

mod.opt              <- final.model.BIC$Opt.model

time.BICEAS          <- proc.time()[3]-st.BICEAS[3]

print("Finally chosen model")
print(mod.opt)

print("Done EAS-BIC")

# MBSP Method (Bai and Ghosh 2018)
st.MBSP              <- proc.time()
final.model.MBSP     <- MBSP::mbsp.tpbn(X = t(X), Y = t(Y))
time.MBSP            <- proc.time()[3]-st.MBSP[3]

print("Done MBSP")

# Sparse Reduced Rank Regression (SRRR) (Chen and Huang (2012))
st.SRRR              <- proc.time()
final.model.SRRR     <- rrpack::srrr(Y = t(Y), X = t(X), nrank = q, method = "adglasso")
time.SRRR            <- proc.time()[3]-st.SRRR[3]

print("Done SRRR")

# Sparse PLS (Chun and Keles 2010)
st.SPLS              <- proc.time()
cv.SPLS              <- spls::cv.spls(x = t(X), y = t(Y), scale.x = F, scale.y = F, eta= seq(0.001,0.99,length.out=20), K = q)
final.model.SPLS     <- spls::spls(x = t(X), y = t(Y), scale.x = F, scale.y = F, eta = cv.SPLS$eta.opt, K = q)
time.SPLS            <- proc.time()[3]-st.SPLS[3]

print("Done SPLS")

# Multivariate Lasso (MLASSO)
st.MLASSO            <- proc.time()
final.model.MLASSO   <- unname(InitialValue_MLasso(Y.m=t(Y), X.m=t(X)))
time.MLASSO          <- proc.time()[3]-st.MLASSO[3]

print("Done MLasso")

# Multivariate Sparse Group Lasso (MSGLASSO) by Li et al. (2015)
st.MSGLASSO           <- proc.time()
final.model.MSGLASSO  <- quiet(unname(InitialValue_MGLasso(Y.m=t(Y), X.m=t(X), lamG=seq(0.01,4,0.05), fold=10)))
time.MSGLASSO         <- proc.time()[3]-st.MSGLASSO[3]

print("Done MSGLasso")

# Dynamic Posterior Exploration with the multivariate spike-and-slab LASSO (mSSL) (Deshpande et al. 2018)
# st.mSSL              <- proc.time()
# final.model.mSSL     <- mSSL::mSSL_dcpe(X = t(X), Y = t(Y))
# time.mSSL            <- proc.time()[3]-st.mSSL[3]

# print("Done DCPE")

# Mutivariate Bayesian Group Lasso with Spike and Slab prior (Liquet et al (2017))
st.MBGLSS            <- proc.time()
final.model.MBGLSS   <- MBSGS::MBGLSS(Y = t(Y), X = t(X), niter=8000, burnin=5000, group_size=rep(1,p), num_update = 10, niter.update = 5, verbose=F)
time.MBGLSS          <- proc.time()[3]-st.MBGLSS[3]

print("Done MBGLSS")

# Out of sample data generation
X.out                <- t(scale(MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = Gamma.mat), center = T, scale = F))
E.out                <- t(MASS::mvrnorm(n = n, mu = rep(0,q), Sigma = Sigma.mat))
Y.out                <- true.B %*% X.out + E.out


# epsilon chosen through OLS based cross-validation and using the B for optimal model
res.CV.OLS_OLS       <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=final.model.CV[[1]]$Opt.B, 
                                    B.true=true.B, model.hat=final.model.CV[[1]]$Opt.model, model.true = true.model,
                                    estModelFreq = final.model.CV[[1]]$ModelFreq)

# epsilon chosen through OLS based cross-validation and using the B by Bayesian model averaging
res.CV.OLS_BMA       <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=final.model.CV[[1]]$Avg.B, 
                                    B.true=true.B, model.hat=final.model.CV[[1]]$Opt.model, model.true = true.model,
                                    estModelFreq = final.model.CV[[1]]$ModelFreq)

# epsilon chosen through Bayesian model averaging based cross-validation and using the B for optimal model
res.CV.BMA_OLS       <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=final.model.CV[[2]]$Opt.B, 
                                    B.true=true.B, model.hat=final.model.CV[[2]]$Opt.model, model.true = true.model,
                                    estModelFreq = final.model.CV[[2]]$ModelFreq)

# epsilon chosen through Bayesian model averaging based cross-validation and using the B by Bayesian model averaging
res.CV.BMA_BMA       <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=final.model.CV[[2]]$Avg.B, 
                                    B.true=true.B, model.hat=final.model.CV[[2]]$Opt.model, model.true = true.model,
                                    estModelFreq = final.model.CV[[2]]$ModelFreq)

# The best model is obtained by BIC for the optimal model and using the B for optimal model
res.BIC.OLS_OLS      <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=final.model.BIC$Opt.B, 
                                    B.true=true.B, model.hat=final.model.BIC$Opt.model, model.true = true.model,
                                    estModelFreq = final.model.BIC$ModelFreq)

# The best model is obtained by BIC for the optimal model and using the B by Bayesian model averaging
res.BIC.OLS_BMA      <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=final.model.BIC$Avg.B, 
                                    B.true=true.B, model.hat=final.model.BIC$Opt.model, model.true = true.model,
                                    estModelFreq = final.model.BIC$ModelFreq)

# The best model is obtained by MBSP developed by Bai and Ghosh (2018)
res.MBSP             <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=t(final.model.MBSP$B.est), 
                                    B.true=true.B, model.hat=final.model.MBSP$active.predictors, model.true = true.model)

# The best model is obtained by Sparse Reduced Rank Regression (Lisha Chen & Jianhua Huang (2012))
B.hat.SRRR           <- t(coef(final.model.SRRR))
model.hat.SRRR       <- which(colMeans(B.hat.SRRR^2) != 0)
res.SRRR             <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.SRRR, 
                                    B.true=true.B, model.hat=model.hat.SRRR, model.true = true.model)

# The best model is obtained by Sparse Penalized Least Squares (Hyonho Chun & Sunduz Keles (2010))
B.hat.SPLS           <- t(coef(final.model.SPLS))
model.hat.SPLS       <- which(colMeans(B.hat.SPLS^2) != 0)
res.SPLS             <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.SPLS, 
                                    B.true=true.B, model.hat=model.hat.SPLS, model.true = true.model)

# The best model is obtained by Multivariate LASSO
B.hat.MLASSO         <- t(final.model.MLASSO)
model.hat.MLASSO     <- which(colMeans(B.hat.MLASSO^2) != 0)
res.MLASSO           <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.MLASSO, 
                                    B.true=true.B, model.hat=model.hat.MLASSO, model.true = true.model)

# The best model is obtained by Multivariate Sparse Group LASSO
B.hat.MSGLASSO       <- t(final.model.MSGLASSO)
model.hat.MSGLASSO   <- which(colMeans(B.hat.MSGLASSO^2) != 0)
res.MSGLASSO         <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.MSGLASSO, 
                                    B.true=true.B, model.hat=model.hat.MSGLASSO, model.true = true.model)

# The best model is obtained by Dynamic Posterior Exploration with the multivariate spike-and-slab LASSO
# B.hat.mSSL           <- t(final.model.mSSL$B)
# model.hat.mSSL       <- which(colMeans(B.hat.mSSL^2) != 0)
# res.mSSL             <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.mSSL, 
#                                     B.true=true.B, model.hat=model.hat.mSSL, model.true = true.model)


# The best model is obtained by Dynamic Posterior Exploration with the multivariate spike-and-slab LASSO
B.hat.MBGLSS         <- t(final.model.MBGLSS$pos_median)
model.hat.MBGLSS     <- which(colMeans(B.hat.MBGLSS^2) != 0)
# To get the average model probabilities
Models.MBGLSS        <- character(length=ncol(final.model.MBGLSS$coef[[1]]))
for (j in 1:ncol(final.model.MBGLSS$coef[[1]])){
  B.hat.MBGLSS.samp  <- sapply(seq_along(final.model.MBGLSS$coef), function(k) final.model.MBGLSS$coef[[k]][,j])
  Models.MBGLSS[j]   <- paste0(which(rowMeans(B.hat.MBGLSS.samp^2) != 0), collapse = " ")
}
res.MBGLSS           <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=t(final.model.MBGLSS$pos_mean), 
                                    B.true=true.B, model.hat=model.hat.MBGLSS, model.true = true.model,
                                    estModelFreq = table(Models.MBGLSS))

res.all              <- cbind(res.CV.OLS_OLS, res.CV.OLS_BMA, res.CV.BMA_OLS, 
                              res.CV.BMA_BMA, res.BIC.OLS_OLS, res.BIC.OLS_BMA, res.MBSP, res.SRRR, 
                              res.SPLS, res.MLASSO, res.MSGLASSO, res.MBGLSS)

res.all              <- rbind(res.all, c(rep(time.CVEAS/4,4), rep(time.BICEAS/2,2), time.MBSP, time.SRRR, 
                                         time.SPLS, time.MLASSO, time.MSGLASSO, time.MBGLSS))


colnames(res.all)    <- c("CV.OLS.OLS", "CV.OLS.BMA", "CV.BMA.OLS", "CV.BMA.BMA", 
                          "BIC.OLS.OLS", "BIC.OLS.BMA", "MBSP", "SRRR", "SPLS", "MLASSO", "MSGLASSO", "MBGLSS")

rownames(res.all)    <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                          "Est.Err.outSample", "Pred.Err.outSample", "FDR", "FNR", "MP", "Prob.TrueModel", 
                          "isTrueModelMAP?", "time (in seconds)")


print(res.all)

assign(paste0("output_f", Exper, "_dsn", dsn), res.all)
save(list=paste0("output_f", Exper, "_dsn", dsn), file=paste0("Results/Experiment", Exper, "/output_f", Exper, "_dsn", dsn, ".Rdata"))
