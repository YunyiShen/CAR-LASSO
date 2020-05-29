library(Matrix)
library(RcppArmadillo)
library(Rcpp)

sourceCpp("./src/sandbox.cpp")

k = 3
n = 4
p = 1


test_data = matrix(rnorm(n * k),k,n)
test_design = matrix(rnorm( n * p ) , n , p)
test_mu = rnorm(k)
test_beta = rnorm(k*p)
test_sigma2 = 1
test_eta2 = runif(k*(k-1))

update_B_helper(test_data,
                test_design,
                test_mu,
                test_beta,
                test_sigma2,
                test_eta2,
                k, n)

