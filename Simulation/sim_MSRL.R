args     <- commandArgs(TRUE)
dsn      <- as.integer(args[1])
n_cores  <- as.integer(args[2])
Exper    <- as.integer(args[3])

print(Exper)
require(Matrix)

cat("Experiment : ", Exper, "\n")
source("setup_BaiGhosh2018.R")
cat("loading dataset :", dsn, "\n")
load(paste0("Datasets/Experiment_", Exper, "/dsn", dsn, ".Rdata"))

source("SourceFunctions/utilities.R")
source("SourceFunctions/Group_MSRL.R")


Y                   <- dat$Y ; 
X                   <- dat$X ; 
true.model          <- dat$true.active ; 
true.B              <- dat$true.B; 


# Out of sample data generation
X.out                <- t(scale(MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = Gamma.mat), center = T, scale = F))
E.out                <- t(MASS::mvrnorm(n = n, mu = rep(0,q), Sigma = Sigma.mat))
Y.out                <- true.B %*% X.out + E.out



# Multivariate square root lasso method (Molstad (2022))
st.MSRL              <- proc.time()
fit.MSRL.cv          <- MSRL.grp.cv(X = t(X), Y = t(Y), nlambda = 50, standardize = FALSE,  
                                    nfolds = 5,  delta = .1, tol = 1e-6, quiet = FALSE, 
                                    inner.quiet = TRUE)
time.MSRL            <- proc.time()[3]-st.MSRL[3]
cat("Done MSRL CV \n")
cat("Time taken by MSRL method : ", time.MSRL, "\n")

final.model.MSRL     <- MSRL.predict(Xnew = t(X.out), fit = fit.MSRL.cv, lambda = fit.MSRL.cv$lam.min)
B.hat.MSRL           <- t(final.model.MSRL$beta)
model.hat.MSRL       <- which(colMeans(B.hat.MSRL^2) != 0)


cat("Done MSRL prediction \n")


# The best model is obtained by MBSP developed by Bai and Ghosh (2018)
res.MSRL             <- get_summary_MSRL(Y.out=Y.out, Y.out.hat = t(final.model.MSRL$pred), B.hat=B.hat.MSRL, 
                                         model.hat=model.hat.MSRL, model.true = true.model)



res.MSRL              <- c(res.MSRL, time.MSRL)


names(res.MSRL)        <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                            "Est.Err.outSample", "Pred.Err.outSample", 
                            "FDR", "FNR", "MP", "Prob.TrueModel", 
                            "isTrueModelMAP?", "time (in seconds)")


print(res.MSRL)

if (!file.exists("Results")){
  dir.create("Results")
}
setwd("Results")
if (!file.exists(paste0('Experiment_', Exper))){
  dir.create(paste0('Experiment_', Exper))
}
setwd(paste0('Experiment_', Exper))
save(res.MSRL, file=paste0("output_MSRL_dsn", dsn, ".Rdata"))