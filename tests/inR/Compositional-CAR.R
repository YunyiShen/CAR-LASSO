library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 8
n = 200
p = 5

sourceCpp("./src/Multinomial-CAR-ALASSO.cpp")
sourceCpp("./src/CAR-ALASSO-hir.cpp")
sourceCpp("./src/CAR-LASSO-hir.cpp")
source("./tests/Formal/Accurancy/Graph_generator.R")

set.seed(42)
Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
#image(Omega)

Sigma <- Graph_raw$Sigma

Design <- matrix(rnorm(n*p,0,1),n,p)
Design <- (Design-mean(Design))/sd(Design)

beta <- as.matrix( rsparsematrix(p,k,0.5))

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




test <- Multinomial_CAR_ALASSO_Cpp(Y,  Design, n_iter = 5000, 
                            n_burn_in = 1000, thin_by = 10, 
                            r_beta = .01+0*beta, delta_beta = 1e-6 + 0 * beta,
                            r_Omega = rep(.01,.5*(k-1)*k),
                            delta_Omega = rep(1e-6,.5*(k-1)*k),
                            lambda_diag = rep(0,k),
                            ns = 100,m = 10, emax = 64,
                            progress = T)

test <- CAR_ALASSO_hir_Cpp(Y,  Design, link = 2,n_iter = 5000, 
                            n_burn_in = 1000, thin_by = 10, 
                            r_beta = .01+0*beta, delta_beta = 1e-6 + 0 * beta,
                            r_Omega = rep(.01,.5*(k-1)*k),
                            delta_Omega = rep(1e-6,.5*(k-1)*k),
                            lambda_diag = rep(0,k),
                            ns = 100,m = 10, emax = 64,
                            progress = T)


test <- CAR_LASSO_hir_Cpp(Y,  Design, link = 2,n_iter = 5000, 
                            n_burn_in = 1000, thin_by = 10, 
                            r_beta = 1, delta_beta = .01,
                            r_Omega = 1,
                            delta_Omega = .01,
                            ns = 100,m = 10, emax = 64,
                            progress = T)

A_Graph <- 0 * Omega
A_Graph[upper.tri(A_Graph,T)] <- apply(test$Omega,2,mean,na.rm = T)
A_Graph <- A_Graph+t(A_Graph)
diag(A_Graph) <- 0.5 * diag(A_Graph)
