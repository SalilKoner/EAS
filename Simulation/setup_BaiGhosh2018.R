# setwd("/Users/salilkoner/Dropbox/EAS_functional/Simulation_Study/SimulationBaiGhosh2018/Sim_For_Revision_EJS")

# devtools::install_github("https://github.com/ajmolstad/MSRL")
# R package MSRL (multivariate square root lasso

if (Exper == 1){
  # Experiment 1  (p < n : sparse model)
  n <- 60 ; p <- 30 ; q <- 3 ; m0 <- 5;
} else if (Exper == 2){
  # Experiment 2  (p < n : dense model)
  n <- 80 ; p <- 60 ; q <- 6 ; m0 <- 40;
} else if (Exper == 3){
  # Experiment 3  (p > n : sparse model)
  n <- 50 ; p <- 200 ; q <- 5 ; m0 <- 20;
} else if (Exper == 4){
  # Experiment 4  (p > n : dense model)
  n <- 60 ; p <- 100 ; q <- 6 ; m0 <- 40;
} else if (Exper == 5){
  # Experiment 5  (p >> n : ultra-sparse model)
  n <- 100 ; p <- 500 ; q <- 3 ; m0 <- 10;
} else if (Exper == 6){
  # Experiment 6  (p >> n : sparse model)
  n <- 150 ; p <- 1000 ; q <- 4 ; m0 <- 50;
} else if (Exper == 7){
  n <- 150; p <- 1000; q <- 60;  m0 <- 50 ; # AE comment: Exploring for large q in high dimensional case
} else if (Exper == 8){
  n <- 150; p <- 1000; q <- 60;  m0 <- 50 # q <- 40; m0 <- 50
} else if (Exper == 9){
  n <- 150; p <- 1000; q <- 60;  m0 <- 50
}


if (Exper %in% 1:7){
  # Covariance of the Design matrix
  Gamma.mat              <- 0.5^{abs(outer(1:p, 1:p,"-"))}
  # Covariance of Noise
  sig                    <- sqrt(2)
  Sigma.mat              <- {0.5^{abs(outer(1:q, 1:q,"-"))}}*{sig^2}
} else if (Exper == 8) {
  # Covariance of the Design matrix
  Gamma.mat              <- toeplitz(c(1, rep(0.5, p-1)))
  # Covariance of Noise
  sig                    <- sqrt(2)
  Sigma.mat              <- {0.5^{abs(outer(1:q, 1:q,"-"))}}*{sig^2}
} else if (Exper == 9) {
  # Covariance of the Design matrix
  Gamma.mat              <- toeplitz(c(1, rep(0.5, p-1)))
  # Covariance of Noise
  sig                    <- sqrt(2)
  Sigma.mat              <- {toeplitz(c(1, rep(0.5, q-1)))}*{sig^2}
}


generate_data          <- function(isim, n, q, p, m0, Sig.X, Sig.E){
  require(MASS)
  X.gen                <- mvrnorm(n = n, mu = rep(0,p), Sigma = Sig.X)
  X                    <- t(scale(X.gen, center=TRUE, scale=FALSE))
  active.index         <- sort(sample(1:p, m0, replace = FALSE))
  B.act                <- runif(m0*q,-5,4)
  B.act                <- matrix(sapply(B.act, function(x) {if (x <= -0.5) return(x) else return(x+1)}),m0,q)
  B                    <- matrix(0, q, p)
  # B[, active.index]    <- matrix(runif(q*m0, 0.5, 5)*{2*rbinom(n = q*m0, size = 1, prob = 0.5)-1}, q, m0)
  B[, active.index]    <- t(B.act)
  E                    <- t(mvrnorm(n = n, mu = rep(0,q), Sigma = Sig.E))
  Y                    <- B %*% X + E
  return(list("Y"=Y, "X"=X, "true.B"=B, "true.active"=active.index))
}

# dat <- generate_data(isim=1, n=n, q=q, p=p, m0=m0, Sig.X=Gamma.mat, Sig.E=Sigma.mat)

