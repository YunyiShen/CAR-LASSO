library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 5
n = 7000
p = 1

sourceCpp("./src/Probit-CAR-ALASSO.cpp")
sourceCpp("./src/Probit-CAR-LASSO.cpp")
sourceCpp("./src/Probit-SRG-LASSO.cpp")
source("./tests/Formal/Accurancy/Graph_generator.R")

set.seed(12345)
#B <- rsparsematrix(k,k,0.2)
#omega <- diag(rgamma(k,1,.1))
#I <- diag(rep(1,k))
#Omega <- t(I-B) %*% omega %*% (I-B)
#diag(Omega) <- diag(Omega) + k
#Omega <- as.matrix(Omega)
Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
#image(Omega)

Sigma <- Graph_raw$Sigma

Design <- matrix(rnorm(n*p,0,1),n,p)
Design <- (Design-mean(Design))/sd(Design)
colnames(Design) <- paste0("x",1:p)


beta <- matrix(rnorm(p*k,0,1),p,k)
#beta[sample(p*k,floor(0.3*p*k))] = 0

mu <- rnorm(k,0,1)
#mu

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)
Y <- Z

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+mu),Sigma)
  Y[i,] <- 1 * ((rnorm(k,Z[i,],1))>0)
}




test <- Probit_CAR_ALASSO_Cpp(Y,  Design, n_iter = 500000, 
                            n_burn_in = 50000, thin_by = 500, 
                            r_beta = 1+0*beta, delta_beta = .01 + 0 * beta,
                            r_Omega = rep(1,.5*(k-1)*k),
                            delta_Omega = rep(.01,.5*(k-1)*k),
                            progress = T)

Graph <- 0 * Omega
Graph[upper.tri(Graph,T)] <- apply(test$Omega,2,mean)
Graph <- Graph+t(Graph)
diag(Graph) <- 0.5 * diag(Graph)

test1 <- Probit_CAR_LASSO_Cpp(Y,  Design, n_iter = 10000, 
                            n_burn_in = 2000, thin_by = 10, 
                            r_beta = 1, delta_beta = .01,
                            r_Omega = 1,
                            delta_Omega = .01,
                            progress = T)

test1 <- Proit_SRG_LASSO_Cpp(Y,  Design, n_iter = 20000, 
                            n_burn_in = 1000, thin_by = 10, 
                            r_beta = 1, delta_beta = .01 ,
                            r_Omega = 1,
                            delta_Omega = .01,
                            progress = T)

Graph1 <- 0 * Omega
Graph1[upper.tri(Graph1,T)] <- apply(test1$Omega,2,mean)
Graph1 <- Graph+t(Graph1)
diag(Graph1) <- 0.5 * diag(Graph1)