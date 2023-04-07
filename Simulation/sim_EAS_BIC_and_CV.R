args     <- commandArgs(TRUE)
dsn      <- as.integer(args[1])
n_cores  <- as.integer(args[2])
Exper    <- as.integer(args[3])
eps.st   <- as.numeric(args[4])
eps.ed   <- as.numeric(args[5])
eps.ln   <- as.integer(args[6])
fold.num <- as.integer(args[7])

cat("Experiment : ", Exper, "\n")
source("setup_BaiGhosh2018.R")
cat("loading dataset :", dsn, "\n")
load(paste0("Datasets/Experiment_", Exper, "/dsn", dsn, ".Rdata"))


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

cat("dimension of Y :", paste(dim(Y), collapse = " by "), "\n")
cat("dimension of X :", paste(dim(X), collapse = " by "), "\n")
cat("dimension of B :", paste(dim(true.B), collapse = " by "), "\n")
cat("length of true.model :", length(true.model), "\n")
cat("Number of grid for epsilon search : ", eps.ln, "\n")
cat("Grid of epsilon search :", eps.grid, "\n") 
cat("Number of folds", fold.num, "\n")

st.CVEAS            <- proc.time()
# EAS method with epsilon chosen through 10-fold cross validation
cvfold              <- admissible_subsets_CV(Y=Y, X=X, N = 100, steps = 500, burnin = 200, epsilonGrid=eps.grid, 
                                             nFold = fold.num, Efficient=TRUE, useMCsample4h=FALSE, addoffset=FALSE, 
                                             parallel=TRUE, n_cores = n_cores, trueModel=true.model, MLasso=TRUE)

# Optimal chosen epsilon
cat("Optimal value of Epsilon chosen by CV method", cvfold$optimal.eps, "\n")

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
cat("Done EAS-CV \n")
cat("time taken by CV method", time.CVEAS, "\n")
cat("Finally chosen model via CV : ", final.model.CV[[1]]$Opt.model, "\n")
cat("True model is : ", true.model)
cat("Number of covariates common in chosen and the true model", 
    length(intersect(final.model.CV[[1]]$Opt.model, true.model)), "\n")

# EAS method with epsilon chosen through BIC
library(foreach); library(doSNOW) ; library(doParallel);

st.BICEAS           <- proc.time()
cl                  <- makeSOCKcluster(n_cores)
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed: %d; tag: %d\n', nfin, tag)) }
opts                <- list(progress=progress)
packages_req        <- c("matrixsampling", "tidyverse", "Rfast", "MSGLasso", "glmnet", "MASS")
BIC.models          <- foreach(epsi=eps.grid, .options.snow=opts, .packages =packages_req) %dopar% {
  admissible_subsets_Efficient(Y=Y, X=X, N=100, steps=2e3, burnin=1e3, PropWeights=NULL,
                               epsilon=epsi, true.Sigma=NULL, addoffset=FALSE,
                               useMCsample4h=FALSE, MLasso=TRUE) }

stopCluster(cl)

BIC.est              <- t(sapply(seq_along(BIC.models), function(i) c(BIC.models[[i]]$Opt.BIC, 
                                                                      BIC.models[[i]]$Avg.BIC, 
                                                                      length(BIC.models[[i]]$Opt.model) )) )
# print(BIC.est)
min.BIC              <- min(BIC.est[,1])
min.BIC.index        <- which(BIC.est[,1] == min.BIC)
eps.chosen           <- eps.grid[min.BIC.index]
cat("Optimal epsilon chosen by BIC method : ", eps.chosen, "\n")

final.model.BIC      <- BIC.models[[sample(min.BIC.index, 1)]]
mod.opt              <- final.model.BIC$Opt.model
time.BICEAS          <- proc.time()[3]-st.BICEAS[3]

cat("Done EAS-BIC\n")
cat("Time taken by BIC method : ", time.BICEAS, "\n")
cat("Finally chosen model via BIC : ", mod.opt, "\n")
cat("True model is : ", true.model, "\n")
cat("Number of covariates common in chosen by BIC and the true model", 
    length(intersect(mod.opt, true.model)), "\n")


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


res.BIC_CV              <- cbind(res.CV.OLS_OLS, res.CV.OLS_BMA, res.CV.BMA_OLS, 
                                 res.CV.BMA_BMA, res.BIC.OLS_OLS, res.BIC.OLS_BMA)

res.BIC_CV              <- rbind(res.BIC_CV, c(rep(time.CVEAS/4,4), rep(time.BICEAS/2,2)))


colnames(res.BIC_CV)    <- c("CV.OLS.OLS", "CV.OLS.BMA", "CV.BMA.OLS", "CV.BMA.BMA", 
                             "BIC.OLS.OLS", "BIC.OLS.BMA")

rownames(res.BIC_CV)    <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                             "Est.Err.outSample", "Pred.Err.outSample", "FDR", "FNR", "MP", "Prob.TrueModel", 
                             "isTrueModelMAP?", "time (in seconds)")


print(res.BIC_CV)

if (!file.exists("Results")){
  dir.create("Results")
}
setwd("Results")
if (!file.exists(paste0('Experiment_', Exper))){
  dir.create(paste0('Experiment_', Exper))
}
setwd(paste0('Experiment_', Exper))
save(res.BIC_CV, file=paste0("output_EAS_BIC_and_CV_dsn", dsn, ".Rdata"))
