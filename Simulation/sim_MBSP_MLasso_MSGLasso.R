args     <- commandArgs(TRUE)
dsn      <- as.integer(args[1])
n_cores  <- as.integer(args[2])
Exper    <- as.integer(args[3])



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


# MBSP Method (Bai and Ghosh 2018)
st.MBSP              <- proc.time()
final.model.MBSP     <- MBSP::MBSP(X = t(X), Y = t(Y), max_steps=6000, burnin=1000, save_samples = FALSE)
time.MBSP            <- proc.time()[3]-st.MBSP[3]

cat("Done MBSP \n")
cat("Time taken by MBSP method : ", time.MBSP, "\n")

# Multivariate Lasso (MLASSO)
st.MLASSO            <- proc.time()
final.model.MLASSO   <- unname(InitialValue_MLasso(Y.m=t(Y), X.m=t(X)))
time.MLASSO          <- proc.time()[3]-st.MLASSO[3]

cat("Done MLasso \n")
cat("Time taken by MLasso method : ", time.MLASSO, "\n")

# Multivariate Sparse Group Lasso (MSGLASSO) by Li et al. (2015)
st.MSGLASSO           <- proc.time()
final.model.MSGLASSO  <- quiet(unname(InitialValue_MGLasso(Y.m=t(Y), X.m=t(X), lamG=seq(0.01,4,0.05), fold=10)))
time.MSGLASSO         <- proc.time()[3]-st.MSGLASSO[3]

cat("Done MSGLasso \n")
cat("Time taken by MSGLasso method : ", time.MSGLASSO, "\n")


# Out of sample data generation
X.out                <- t(scale(MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = Gamma.mat), center = T, scale = F))
E.out                <- t(MASS::mvrnorm(n = n, mu = rep(0,q), Sigma = Sigma.mat))
Y.out                <- true.B %*% X.out + E.out


cat("Extracting summary for MBSP \n")

# The best model is obtained by MBSP developed by Bai and Ghosh (2018)
res.MBSP             <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=t(final.model.MBSP$B_est), 
                                    B.true=true.B, model.hat=final.model.MBSP$active_predictors, model.true = true.model)


cat("Extracting summary for MLasso \n")

# The best model is obtained by Multivariate LASSO
B.hat.MLASSO         <- t(final.model.MLASSO)
model.hat.MLASSO     <- which(colMeans(B.hat.MLASSO^2) != 0)
res.MLASSO           <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.MLASSO, 
                                    B.true=true.B, model.hat=model.hat.MLASSO, model.true = true.model)


cat("Extracting summary for MSGLasso \n")

# The best model is obtained by Multivariate Sparse Group LASSO
B.hat.MSGLASSO       <- t(final.model.MSGLASSO)
model.hat.MSGLASSO   <- which(colMeans(B.hat.MSGLASSO^2) != 0)
res.MSGLASSO         <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.MSGLASSO, 
                                    B.true=true.B, model.hat=model.hat.MSGLASSO, model.true = true.model)


res.MBSP_MLasso_MSGLasso  <- cbind(res.MBSP, res.MLASSO, res.MSGLASSO)

res.MBSP_MLasso_MSGLasso  <- rbind(res.MBSP_MLasso_MSGLasso, c(time.MBSP, time.MLASSO, time.MSGLASSO))


colnames(res.MBSP_MLasso_MSGLasso)    <- c("MBSP", "MLASSO", "MSGLASSO")

rownames(res.MBSP_MLasso_MSGLasso)    <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                                           "Est.Err.outSample", "Pred.Err.outSample", "FDR", "FNR", "MP", "Prob.TrueModel", 
                                           "isTrueModelMAP?", "time (in seconds)")


print(res.MBSP_MLasso_MSGLasso)

if (!file.exists("Results")){
  dir.create("Results")
}
setwd("Results")
if (!file.exists(paste0('Experiment_', Exper))){
  dir.create(paste0('Experiment_', Exper))
}
setwd(paste0('Experiment_', Exper))
save(res.MBSP_MLasso_MSGLasso, file=paste0("output_MBSP_MLasso_MSGLasso_dsn", dsn, ".Rdata"))
