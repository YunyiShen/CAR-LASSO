library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

sourceCpp("./src/Probit_helper.cpp")

k = 2
n = 2
p = 1

source("./tests/Formal/Accurancy/Graph_generator.R")

set.seed(12345)

Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
#image(Omega)

Sigma <- Graph_raw$Sigma

mu <- c(2,-2)
mu_curr <- as.matrix(mu)
Design <- matrix(1,2,1)
beta_curr <- matrix(0,p,k)

Z <- matrix(NA,n,k)
Y <- Z

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% mu,Sigma)
  Y[i,] <- 1 * ((rnorm(k,Z[i,],1))>0)
}



nrep <- 1000
res <- matrix(NA,nrep,n*k)
Z_curr <- matrix(0,n,k)



for(i in 1:nrep){
    update_Z_helper_CAR(Z_curr,Y,Design,mu_curr,beta_curr,Omega,k,p,n)
    res[i,] <- c(Z_curr)

}
