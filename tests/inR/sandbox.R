library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())
sourceCpp("./src/sandbox.cpp")

k = 5
n = 20
p = 2

test_sigma <- matrix(rnorm(k^2),k,k)
test_sigma <- t(test_sigma) %*% test_sigma
test_data = MASS::mvrnorm(n,c(0,0,0,0,0),test_sigma)

test_design <- matrix(rnorm(n*p),n,p)

test_mu <- rnorm(k)
test_tau2 <- runif(k*p)

test_beta = matrix(rnorm(k*p),p,k)

test_lambda_Omega <- runif(1)

#update_beta_helper(test_data,test_design,test_mu,test_tau2,test_sigma,k,p,n)
#update_mu_helper(test_data,test_design,test_beta,test_sigma,k,p,n)
update_tau2_helper(test_beta,
                    test_lambda_Omega,
                    test_sigma,
                    k,p,n)
