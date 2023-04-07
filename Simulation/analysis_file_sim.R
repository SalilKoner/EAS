
args     <- commandArgs(TRUE)
Exper    <- as.integer(args[1])
n_rep    <- as.integer(args[2])

print(paste0("Working for experiment: ", Exper))

sim.result <- array(NA,dim=c(11,6,n_rep))

for (dsn in 1:n_rep){
  filename         <-  paste0("Results_EAS/Experiment_", Exper, "/output_EAS_BIC_and_CV_dsn", dsn, ".Rdata")
  if (file.exists(filename)){
    load(filename)
    sim.result[, , dsn] <- res.BIC_CV
  } else{
    cat("file", dsn, "does not exist \n")
  }
}

Summary.mean   <- round(apply(sim.result, c(1,2), mean, na.rm=T),4)
Summary.median <- round(apply(sim.result, c(1,2), median, na.rm=T),4)

colnames(Summary.mean)    <- colnames(Summary.median) <- c("CV.OLS.OLS", "CV.OLS.BMA", "CV.BMA.OLS", "CV.BMA.BMA", 
                                                           "BIC.OLS.OLS", "BIC.OLS.BMA") 

rownames(Summary.mean)    <- rownames(Summary.median)   <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                                                             "Est.Err.outSample", "Pred.Err.outSample", 
                                                             "FDR", "FNR", "MP", "Prob.TrueModel", "isTrueModelMAP?",
                                                             "time (in seconds)")

#save(list=c("Summary.mean", "Summary.median"), file=paste0("Results/Experiment_", Exper, "/Summary.Rdata"))


# print(t(Summary.mean)[,c(1,3,5,6:10)])
# print(t(Summary.median)[,c(1,3,5,6:10)])

print(cbind(t(Summary.median)[,c(5, 11), drop=F], t(Summary.mean)[,c(6:10), drop=F]))




# MBGLSS + MLASSO + MSGLasso

args     <- commandArgs(TRUE)
Exper    <- as.integer(args[1])

Exper <- 9

print(paste0("Working for experiment: ", Exper))

sim.result <- array(NA,dim=c(11,3,n_rep))

for (dsn in 1:n_rep){
  filename         <-  paste0("Results_other_methods/Experiment_", Exper, "/MBSP_MLasso_MSGLasso/output_MBSP_MLasso_MSGLasso_dsn", dsn, ".Rdata")
  if (file.exists(filename)){
    load(filename)
    sim.result[, , dsn] <- res.MBSP_MLasso_MSGLasso
  } else{
    cat("file", dsn, "does not exist \n")
  }
}

Summary.mean   <- round(apply(sim.result, c(1,2), mean, na.rm=T),4)
Summary.median <- round(apply(sim.result, c(1,2), median, na.rm=T),4)

colnames(Summary.mean)    <- colnames(Summary.median) <- c("MBSP", "MLasso", "MSGLasso")

rownames(Summary.mean)    <- rownames(Summary.median)   <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                                                             "Est.Err.outSample", "Pred.Err.outSample", 
                                                             "FDR", "FNR", "MP", "Prob.TrueModel", "isTrueModelMAP?",
                                                             "time (in seconds)")

#save(list=c("Summary.mean", "Summary.median"), file=paste0("Results/Experiment_", Exper, "/Summary.Rdata"))


print(t(Summary.mean)[,c(1,3,5,6:10)])
print(t(Summary.median)[,c(1,3,5,6:10)])

print(cbind(t(Summary.median)[,c(5,11), drop=F], t(Summary.mean)[,c(6:10), drop=F]))


# SPLS and SRRR
args     <- commandArgs(TRUE)
Exper    <- as.integer(args[1])

Exper    <- 9

print(paste0("Working for experiment: ", Exper))

sim.result <- array(NA,dim=c(11,2,n_rep))

for (dsn in 1:n_rep){
  filename         <-  paste0("Results_other_methods/Experiment_", 
                              Exper, "/SPLS_SRRR/output_SPLS_SRRR_dsn", dsn, ".Rdata")
  if (file.exists(filename)){
    load(filename)
    sim.result[, , dsn] <- res.SPLS_SRRR
  } else{
    cat("file", dsn, "does not exist \n")
  }
}

Summary.mean   <- round(apply(sim.result, c(1,2), mean, na.rm=T),4)
Summary.median <- round(apply(sim.result, c(1,2), median, na.rm=T),4)

colnames(Summary.mean)    <- colnames(Summary.median) <- c("SPLS", "SRRR")

rownames(Summary.mean)    <- rownames(Summary.median)   <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                                                             "Est.Err.outSample", "Pred.Err.outSample", 
                                                             "FDR", "FNR", "MP", "Prob.TrueModel", "isTrueModelMAP?",
                                                             "time (in seconds)")

#save(list=c("Summary.mean", "Summary.median"), file=paste0("Results/Experiment_", Exper, "/Summary.Rdata"))


# print(t(Summary.mean)[,c(1,3,5,6:10)])
# print(t(Summary.median)[,c(1,3,5,6:10)])

print(cbind(t(Summary.median)[,c(5,11), drop=F], t(Summary.mean)[,c(6:10, drop=F)]))



# MSRL

args     <- commandArgs(TRUE)
Exper    <- as.integer(args[1])

Exper    <- 9

print(paste0("Working for experiment: ", Exper))

sim.result <- array(NA,dim=c(11,1,n_rep))

for (dsn in 1:n_rep){
  filename         <-  paste0("Results_other_methods/Experiment_", 
                              Exper,  "/MSRL/output_MSRL_dsn", dsn, ".Rdata")
  if (file.exists(filename)){
    load(filename)
    sim.result[, , dsn] <- res.MSRL
  } else{
    cat("file", dsn, "does not exist \n")
  }
}

Summary.mean   <- round(apply(sim.result, c(1,2), mean, na.rm=T),4)
Summary.median <- round(apply(sim.result, c(1,2), median, na.rm=T),4)

colnames(Summary.mean)    <- colnames(Summary.median) <- c("MSRL")

rownames(Summary.mean)    <- rownames(Summary.median)   <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                                                             "Est.Err.outSample", "Pred.Err.outSample", 
                                                             "FDR", "FNR", "MP", "Prob.TrueModel", "isTrueModelMAP?",
                                                             "time (in seconds)")

#save(list=c("Summary.mean", "Summary.median"), file=paste0("Results/Experiment_", Exper, "/Summary.Rdata"))

# print(t(Summary.mean)[,c(1,3,5,6:10)])
# print(t(Summary.median)[,c(1,3,5,6:10)])

print(cbind(t(Summary.median)[,c(5,11), drop=F], t(Summary.mean)[,c(6:10), drop=F]))

