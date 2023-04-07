args     <- commandArgs(TRUE)
dsn      <- as.integer(args[1])
n_cores  <- as.integer(args[2])
Exper    <- as.integer(args[3])



cat("Experiment : ", Exper, "\n")
source("setup_BaiGhosh2018.R")
# load data
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


# Sparse Reduced Rank Regression (SRRR) (Chen and Huang (2012))
st.SRRR              <- proc.time()
final.model.SRRR     <- rrpack::srrr(Y = t(Y), X = t(X), nrank = q)
time.SRRR            <- proc.time()[3]-st.SRRR[3]

cat("Done SRRR \n")
cat("Time taken by SRRR method : ", time.SRRR, "\n")

# Sparse PLS (Chun and Keles 2010)
st.SPLS              <- proc.time()
cv.SPLS              <- spls::cv.spls(x = t(X), y = t(Y), scale.x = F, scale.y = F, 
                                      eta= seq(0.01,0.99,length.out=10), K = q)
final.model.SPLS     <- spls::spls(x = t(X), y = t(Y), scale.x = F, scale.y = F,
                                   eta = cv.SPLS$eta.opt, K = cv.SPLS$K.opt)
time.SPLS            <- proc.time()[3]-st.SPLS[3]

cat("Done SPLS \n")
cat("Time taken by SPLS method : ", time.SPLS, "\n")


# Out of sample data generation
X.out                <- t(scale(MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = Gamma.mat), center = T, scale = F))
E.out                <- t(MASS::mvrnorm(n = n, mu = rep(0,q), Sigma = Sigma.mat))
Y.out                <- true.B %*% X.out + E.out


cat("Extracting summary for SRRR \n")

# The best model is obtained by Sparse Reduced Rank Regression (Lisha Chen & Jianhua Huang (2012))
B.hat.SRRR           <- t(coef(final.model.SRRR))
model.hat.SRRR       <- which(colMeans(B.hat.SRRR^2) != 0)
res.SRRR             <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.SRRR, 
                                    B.true=true.B, model.hat=model.hat.SRRR, model.true = true.model)


cat("Extracting summary for SPLS \n")

# The best model is obtained by Sparse Penalized Least Squares (Hyonho Chun & Sunduz Keles (2010))
B.hat.SPLS           <- t(coef(final.model.SPLS))
model.hat.SPLS       <- which(colMeans(B.hat.SPLS^2) != 0)
res.SPLS             <- get_summary(Y.in=Y, X.in=X, Y.out=Y.out, X.out=X.out, B.hat=B.hat.SPLS, 
                                    B.true=true.B, model.hat=model.hat.SPLS, model.true = true.model)


res.SPLS_SRRR       <- cbind(res.SPLS, res.SRRR)

res.SPLS_SRRR       <- rbind(res.SPLS_SRRR, c(time.SPLS, time.SRRR))


colnames(res.SPLS_SRRR)    <- c("SPLS", "SRRR")

rownames(res.SPLS_SRRR)    <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                                     "Est.Err.outSample", "Pred.Err.outSample", 
                                     "FDR", "FNR", "MP", "Prob.TrueModel", 
                                     "isTrueModelMAP?", "time (in seconds)")


print(res.SPLS_SRRR)

if (!file.exists("Results")){
  dir.create("Results")
}
setwd("Results")
if (!file.exists(paste0('Experiment_', Exper))){
  dir.create(paste0('Experiment_', Exper))
}
setwd(paste0('Experiment_', Exper))
save(res.SPLS_SRRR, file=paste0("output_SPLS_SRRR_dsn", dsn, ".Rdata"))
