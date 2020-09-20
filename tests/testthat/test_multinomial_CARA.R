test_that("Multinomial-CAR-ALASSO.cpp",{
    library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 5
n = 8000
p = 1

source("./tests/Formal/Accurancy/Graph_generator.R")

set.seed(12345)
Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
#image(Omega)

Sigma <- Graph_raw$Sigma

Design <- matrix(rnorm(n*p,0,1),n,p)
Design <- (Design-mean(Design))/sd(Design)

beta <- matrix(rnorm(p*k,0,1),p,k)

mu <- rnorm(k,0,1)

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)
Y <- matrix(NA,n,k+1)

N = 5000

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+mu),Sigma)
  p_temp <- c(exp(Z[i,]) , 1)
  Y[i,] <- rmultinom(1, N, p_temp/sum(p_temp))
}




test <- Multinomial_CAR_ALASSO_Cpp(Y,  Design, n_iter = 50000, 
                            n_burn_in = 5000, thin_by = 50, 
                            r_beta = 1+0*beta, delta_beta = .01 + 0 * beta,
                            r_Omega = rep(1,.5*(k-1)*k),
                            delta_Omega = rep(.01,.5*(k-1)*k),
                            ns = 5000,m = 100, emax = 640,
                            progress = T)

})