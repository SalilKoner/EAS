
args     <- commandArgs(TRUE)
Exper    <- as.integer(args[1])
n_rep    <- as.integer(args[2])

print(paste0("Working for experiment: ", Exper))

sim.result <- array(NA,dim=c(11,ifelse(Exper < 7,9,8),n_rep))

for (dsn in 1:n_rep){
  filename         <-  paste0("Results/Experiment_", Exper, 
                              "/output_EAS_BIC_and_CV_dsn", dsn, ".Rdata")
  if (file.exists(filename)){
    load(filename)
    sim.result[, 1:2, dsn] <- res.BIC_CV
  } else{
    cat("file", dsn, "does not exist \n")
  }
  
  filename         <-  paste0("Results/Experiment_", Exper, 
                              "/output_MBSP_MLasso_MSGLasso_dsn", dsn, ".Rdata")
  if (file.exists(filename)){
    load(filename)
    sim.result[, 3:5, dsn] <- res.MBSP_MLasso_MSGLasso
  } else{
    cat("file", dsn, "does not exist \n")
  }
  
  filename         <-  paste0("Results/Experiment_", 
                              Exper, "/output_SPLS_SRRR_dsn", dsn, ".Rdata")
  if (file.exists(filename)){
    load(filename)
    sim.result[, 6:7, dsn] <- res.SPLS_SRRR
  } else{
    cat("file", dsn, "does not exist \n")
  }
  
  filename         <-  paste0("Results/Experiment_", 
                              Exper,  "/output_MSRL_dsn", dsn, ".Rdata")
  if (file.exists(filename)){
    load(filename)
    sim.result[, 8, dsn] <- res.MSRL
  } else{
    cat("file", dsn, "does not exist \n")
  }
  
  if (Exper %in% 1:6){
    filename         <-  paste0("Results/Experiment_", 
                                Exper,  "/output_MBGLSS_dsn", dsn, ".Rdata")
    if (file.exists(filename)){
      load(filename)
      sim.result[, 9, dsn] <- res.MBGLSS
    } else{
      cat("file", dsn, "does not exist \n")
    }
  }

}

Summary.mean   <- round(apply(sim.result, c(1,2), mean, na.rm=T),4)
Summary.median <- round(apply(sim.result, c(1,2), median, na.rm=T),4)

colnames(Summary.mean)    <- colnames(Summary.median) <- c(c("CV.OLS.OLS", "BIC.OLS.OLS",
                                                             "MBSP", "MLasso", "MSGLasso",
                                                             "SPLS", "SRRR",
                                                             "MSRL"), switch(Exper < 7, "MBGLSS")) 

rownames(Summary.mean)    <- rownames(Summary.median)   <- c("Err.B", "Est.Err.inSample", "Pred.Err.inSample",
                                                             "Est.Err.outSample", "Pred.Err.outSample", 
                                                             "FDR", "FNR", "MP", "Prob.TrueModel", "isTrueModelMAP?",
                                                             "time (in seconds)")

cat("Results for Experiment - ", Exper, ": \n")

result  <- cbind(t(Summary.median)[,c(5), drop=F], 
                 t(Summary.mean)[,c(6:10, drop=F)],
                 t(Summary.median)[,c(11), drop=F])

print(result)

save(result, file=paste0("Results/Experiment_", Exper, "/Summary_all", ".Rdata"))
