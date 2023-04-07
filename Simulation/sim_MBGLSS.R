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

# Mutivariate Bayesian Group Lasso with Spike and Slab prior (Liquet et al (2017))
st.MBGLSS            <- proc.time()
final.model.MBGLSS   <- MBSGS::MBGLSS(Y = t(Y), X = t(X), niter=8000, burnin=5000, group_size=rep(1,p), num_update = 10, niter.update = 5, verbose=F)
time.MBGLSS          <- proc.time()[3]-st.MBGLSS[3]

print("Done MBGLSS")

# Out of sample data generation
X.out                <- t(scale(MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = Gamma.mat), center = T, scale = F))
E.out                <- t(MASS::mvrnorm(n = n, mu = rep(0,q), Sigma = Sigma.mat))
Y.out                <- true.B %*% X.out + E.out


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

res.MBGLSS              <- c(res.MBGLSS, time.MBGLSS)


names(res.MBGLSS)        <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                            "Est.Err.outSample", "Pred.Err.outSample", 
                            "FDR", "FNR", "MP", "Prob.TrueModel", 
                            "isTrueModelMAP?", "time (in seconds)")


print(res.MBGLSS)

if (!file.exists("Results")){
  dir.create("Results")
}
setwd("Results")
if (!file.exists(paste0('Experiment_', Exper))){
  dir.create(paste0('Experiment_', Exper))
}
setwd(paste0('Experiment_', Exper))
save(res.MBGLSS, file=paste0("output_MBGLSS_dsn", dsn, ".Rdata"))