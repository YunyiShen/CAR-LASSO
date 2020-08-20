library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 10
n = 1000
p = 2


sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")
sourceCpp("./src/SRG-LASSO.cpp")
source("./R/misc.R")

B <- rsparsematrix(k,k,0.2)
omega <- diag(rgamma(k,3,.1))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
#Omega <- omega
#diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)
image(Omega)

Sigma <- solve(Omega)

Design <- 1.0* (matrix(rnorm(n*p,0,1),n,p))
#Design <- (Design-mean(Design))/sd(Design)
colnames(Design) <- paste0("x",1:p)


beta <- matrix(rnorm(p*k,2,1),p,k)
#beta[sample(p*k,floor(0.3*p*k))] = 0

mu <-  1+rnorm(k)
#mu

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)
Y <- Z

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+mu),Sigma)
  Y[i,] <- 1 * ((rnorm(k,Z[i,],1))>0)
}

par(mfrow = c(1,2))

CAR_test <- CAR_LASSO_Cpp(Z,  Design, n_iter = 25000, 
                          n_burn_in = 5000, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)


CAR_Graph <- 0 * Omega
CAR_Graph[upper.tri(CAR_Graph,T)] <- apply(CAR_test$Omega,2,mean)
CAR_Graph <- CAR_Graph+t(CAR_Graph)
diag(CAR_Graph) <- 0.5 * diag(CAR_Graph)
image((CAR_Graph))
image(Omega)
hist((CAR_Graph-Omega)/Omega)


SRG_test <- SRG_LASSO_Cpp(Z,  Design, n_iter = 25000, 
                          n_burn_in = 5000, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)


SRG_Graph <- 0 * Omega
SRG_Graph[upper.tri(SRG_Graph,T)] <- apply(SRG_test$Omega,2,mean)
SRG_Graph <- SRG_Graph+t(SRG_Graph)
diag(SRG_Graph) <- 0.5 * diag(SRG_Graph)




Glasso <- Graphical_LASSO_Cpp(Z, 25000, 5000, 10, 1, .01, T)

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
SRG_stein_loss <- stein_loss_cpp(SRG_Graph,Omega)
Glasso_stein_loss <- stein_loss_cpp(Glasso_Graph,Omega)

