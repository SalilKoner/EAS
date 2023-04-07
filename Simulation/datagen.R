

args    <- commandArgs(TRUE)
Exper   <- as.integer(args[1])
n_cores <- as.integer(args[2])
n_data  <- as.integer(args[3])

cur.dir <- getwd()

cat("working directory", cur.dir, "\n")
cat("Experiment", Exper, "\n")
source("setup_BaiGhosh2018.R")
cat("n = ", n, "p = ", p, "q = ", q, "m0 = ", m0, "\n")
cat("number of replications : ", n_data, "\n")

if (!file.exists("Datasets")){
  dir.create("Datasets")
}
setwd("Datasets")
if (!file.exists(paste0('Experiment_', Exper))){
  dir.create(paste0('Experiment_', Exper))
}
setwd(paste0('Experiment_', Exper))
for (dsn in 1:n_data){
  dat <- generate_data(isim=dsn, n=n, q=q, p=p, m0=m0, Sig.X=Gamma.mat, Sig.E=Sigma.mat)
  save(dat, file=paste0("dsn", dsn, ".Rdata"))
}