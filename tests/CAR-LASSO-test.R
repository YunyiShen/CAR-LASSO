library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 11
n = 8000
p = 2


sourceCpp("./src/CAR-LASSO.cpp")

B <- rsparsematrix(k,k,0.3)
omega <- diag(rgamma(k,10,.3))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
#diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)
image(Omega)

Sigma <- solve(Omega)

Design <- matrix(rnorm(n*p,0,1),n,p)
Design <- (Design-mean(Design))/sd(Design)
colnames(Design) <- paste0("x",1:p)


beta <- matrix(rnorm(p*k,1,1),p,k)
#beta[sample(p*k,floor(0.3*p*k))] = 0

mu <-  rnorm(k)
#mu

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)
Y <- Z

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+mu),Sigma)
  Y[i,] <- 1 * ((rnorm(k,Z[i,],1))>0)
}


CAR_test <- CAR_LASSO_Cpp(Z,  Design, n_iter = 5000, 
                          n_burn_in = 1000, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)


CAR_Graph <- 0 * Omega
CAR_Graph[upper.tri(CAR_Graph,T)] = apply(CAR_test$Omega,2,median)
CAR_Graph = CAR_Graph+t(CAR_Graph)
diag(CAR_Graph) = 0.5 * diag(CAR_Graph)
image((CAR_Graph))
image(Omega)
hist((CAR_Graph-Omega)/Omega)





