library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 8
n = 7466
p = 1

sourceCpp("./src/SRG-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")

set.seed(42)
Sigma <- matrix(rnorm(k^2,0,3),k,k)
Sigma <- t(Sigma) %*% Sigma

Design <- matrix(rnorm(n*p),n,p)
beta <- matrix(rnorm(1),p,k)
mu <- rnorm(k)

Xbeta <- Design %*% beta

simu <- matrix(NA,n,k)

for( i in 1:n ){
  simu[i,] <- MASS::mvrnorm(1,Xbeta[i,]+mu,Sigma)
}

test <- SRG_LASSO_Cpp(simu,  Design, n_iter = 10000, 
                      n_burn_in = 500, thin_by = 10, 
                      r_beta = 1, delta_beta = .1,
                      r_Omega = 1,delta_Omega = .1,
                      progress = T)

Glasso <- Graphical_LASSO_Cpp(simu,2000,1000,10,1,1,T)
