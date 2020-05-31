library(Matrix)
library(RcppArmadillo)
library(Rcpp)

rm(list = ls())
sourceCpp("./src/sandbox.cpp")

k = 3
n = 4
p = 1


test_data = matrix(rnorm(n * k),k,n)
test_design = matrix(rnorm( n * p ) , n , p)
test_mu = rnorm(k)
test_beta = rnorm(k*p)
test_B = rnorm(k*(k-1))
test_sigma2 = 1
test_tau2 = runif(k*p)
test_eta2 = runif(k*(k-1))
test_lambda2_beta = 1
test_lambda2_B = 1


test_sigma_vec = runif(k)
test_B = matrix(rnorm(k*k),k,k)
diag(test_B) = 0


set.seed(42)
test_SAR = Sample_SAR_Cpp(10000,test_mu,test_sigma_vec,test_B)

cov(t(test_SAR))
rowMeans(test_SAR)
