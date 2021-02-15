library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 10
n = 200
p = 10


sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")
sourceCpp("./src/CAR-ALASSO.cpp")
source("./R/misc.R")

source("./tests/Formal/Accurancy/Graph_generator.R")

set.seed(12345)
Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
#image(Omega)

Sigma <- Graph_raw$Sigma

#B <- rsparsematrix(k,k,0.2)
#omega <- diag(rgamma(k,1,.1))
#I <- diag(rep(1,k))
#Omega <- t(I-B) %*% omega %*% (I-B)
#diag(Omega) <- diag(Omega) + k
#Omega <- as.matrix(Omega)

Design <- 1.0* (matrix(rnorm(n*p,0,1),n,p))
#Design <- (Design-mean(Design))/sd(Design)
colnames(Design) <- paste0("x",1:p)


beta <- matrix(rnorm(p*k,0,5),p,k)
beta[sample(p*k,floor(0.5*p*k))] = 0

mu <-  rnorm(k)
#mu

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+mu),Sigma)
}

par(mfrow = c(1,2))

CAR_test <- CAR_LASSO_Cpp(Z,  Design, n_iter = 10000, 
                          n_burn_in = 5000, thin_by = 10, 
                          r_beta = 1, delta_beta = 1e-2,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)

CAR_Graph <- 0 * Omega
CAR_Graph[upper.tri(CAR_Graph,T)] <- apply(CAR_test$Omega,2,mean)
CAR_Graph <- CAR_Graph+t(CAR_Graph)
diag(CAR_Graph) <- 0.5 * diag(CAR_Graph)
image((CAR_Graph))
image(Omega)
hist((CAR_Graph-Omega)/Omega)


CAR_A_test <- CAR_ALASSO_Cpp(Z,  Design, n_iter = 10000, 
                          n_burn_in = 5000, thin_by = 10, 
                          r_beta = 1e-2+0*beta, delta_beta = 1e-6 + 0 * beta,
                          r_Omega = rep(1e-2,.5*(k-1)*k),
                          delta_Omega = rep(1e-6,.5*(k-1)*k),
                          lambda_diag = rep(0,k),
                          progress = T)

A_Graph <- 0 * Omega
A_Graph[upper.tri(A_Graph,T)] <- apply(CAR_A_test$Omega,2,mean)
A_Graph <- A_Graph+t(A_Graph)
diag(A_Graph) <- 0.5 * diag(A_Graph)


multireg <- CAR_multireg(Z,Design,10000)

multireg_Graph <- 0 * Omega
multireg_Graph[upper.tri(multireg_Graph,T)] <- apply(multireg$Omega,2,mean)
multireg_Graph <- multireg_Graph+t(multireg_Graph)
diag(multireg_Graph) <- 0.5 * diag(multireg_Graph)


Glasso <- Graphical_LASSO_Cpp(Z, 10000, 5000, 10, 1, .01, T)

Glasso_Graph <- 0 * Omega
Glasso_Graph[upper.tri(Glasso_Graph,T)] = apply(Glasso$Omega,2,mean)
Glasso_Graph <- Glasso_Graph+t(Glasso_Graph)
diag(Glasso_Graph) <- 0.5 * diag(Glasso_Graph)
image(Glasso_Graph)


image((abs(Glasso_Graph)>.1*sd(Glasso_Graph[upper.tri(Omega)])))
image((abs(CAR_Graph)>.1*sd(CAR_Graph[upper.tri(Omega)])))

image((abs(Omega)>0*sd(Omega[upper.tri(Omega)]))-
        (abs(CAR_Graph)>.2*sd(CAR_Graph[upper.tri(Omega)])))

image((abs(Omega)>0*sd(Omega[upper.tri(Omega)]))-
        (abs(Glasso_Graph)>.2*sd(Glasso_Graph[upper.tri(Omega)])))


hist((Glasso_Graph-Omega)/Omega)
image(Omega)

mean(((Glasso_Graph-Omega))^2)


image(CAR_Graph-Glasso_Graph)
hist(CAR_Graph-Glasso_Graph)

mean(((CAR_Graph-Omega)^2))


CAR_stein_loss <- stein_loss_cpp(CAR_Graph,Omega)
A_stein_loss <- stein_loss_cpp(A_Graph,Omega)
Glasso_stein_loss <- stein_loss_cpp(Glasso_Graph,Omega)

