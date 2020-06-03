library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

sourceCpp("./src/graphical-LASSO.cpp")

test_sigma <- matrix(rnorm(25),5,5)
test_sigma <- t(test_sigma) %*% test_sigma
test_mat = MASS::mvrnorm(5000,c(0,0,0,0,0),test_sigma)

set.seed(42)
ttt = blockGLasso(test_mat,iterations = 5000)

Omega_post <- Reduce(`+`,ttt$Omegas)/length(ttt$Omegas)

set.seed(42)
tt_cpp = Graphical_LASSO_Cpp(test_mat,5000,1000,100,1,0.1,T)


Omega_postcpp <- colMeans(tt_cpp$Omega)
Omega_postcpp = matrix(Omega_postcpp,5,5)


Omega_post-Omega_postcpp
