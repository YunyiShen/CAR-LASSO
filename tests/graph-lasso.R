library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

sourceCpp("./src/graphical-LASSO.cpp")

test_sigma <- matrix(rnorm(25),5,5)
test_sigma <- t(test_sigma) %*% test_sigma
Omega <- solve(test_sigma)

test_mat = MASS::mvrnorm(2000,c(0,0,0,0,0),test_sigma)

#set.seed(42)
ttt = blockGLasso(test_mat,iterations = 50000)

Omega_post <- Reduce(`+`,ttt$Omegas)/length(ttt$Omegas)
(Omega-Omega_post)/Omega


  set.seed(42)
  test_sigma <- matrix(rnorm(25),5,5)
  test_sigma <- t(test_sigma) * test_sigma
  Omega <- solve(test_sigma)
  
  test_mat = MASS::mvrnorm(8000,rep(0,5),test_sigma)
  BGL_package_res <- BayesianGLasso::blockGLasso(test_mat,iterations = 5000)
  Omega_post_pkg <- Reduce(`+`,BGL_package_res$Omegas)/length(BGL_package_res$Omegas)
  (Omega-Omega_post_pkg)/Omega

#set.seed(42)
tt_cpp = Graphical_LASSO_Cpp(test_mat,50000,1000,100,1,0.1,T)


Omega_posttemp <- colMeans(tt_cpp$Omega)
Omega_postcpp = matrix(0,5,5)
Omega_postcpp[upper.tri(Omega_postcpp,T)] = Omega_posttemp
Omega_postcpp = Omega_postcpp+t(Omega_postcpp)
diag(Omega_postcpp) = diag(Omega_postcpp)/2

(Omega_post-Omega_post_pkg)/Omega_post
