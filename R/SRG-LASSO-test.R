library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 11
n = 8000
p = 1

sourceCpp("./src/Probit-SRG-LASSO.cpp")
sourceCpp("./src/SRG-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")


#set.seed(42)
Omega <- rsparsematrix(k,k,0.5)
Omega <- Omega*t(Omega)
diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)

Sigma <- solve(Omega)

Design <- matrix(rnorm(n*p,0,3),n,p)
Design <- (Design-mean(Design))/sd(Design)
beta <- matrix(rnorm(p*k,1.1,0.05),p,k)
mu <- rnorm(k)
#mu

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)
Y <- Z

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Xbeta[i,]+mu,Sigma)
  Y[i,] <- 1 * ((rnorm(k,Z[i,],1))>0)
}




test <- Proit_SRG_LASSO_Cpp(Y,  Design, n_iter = 2000, 
                            n_burn_in = 1000, thin_by = 10, 
                            r_beta = 1, delta_beta = .1,
                            r_Omega = 1,delta_Omega = .1,
                            progress = T)

#par(mfrow = c(1,2))
Omega_uptri <- Omega[upper.tri(Omega,T)]
diff_probit <- apply(test$Omega,1,function(w,k){(w-k)/k},Omega_uptri)
hist(diff_probit)

SRG_test <- SRG_LASSO_Cpp(Y,  Design, n_iter = 2000, 
                          n_burn_in = 1000, thin_by = 10, 
                          r_beta = 1, delta_beta = .1,
                          r_Omega = 1,delta_Omega = .1,
                          progress = T)

diff_SRG <- apply(SRG_test$Omega,1,function(w,k){(w-k)/k},Omega_uptri)
hist(diff_SRG)

Glasso <- Graphical_LASSO_Cpp(Z,2000,1000,10,1,.1,T)
diff_Glasso <- apply(Glasso$Omega,1,function(w,k){(w-k)/k},c(Omega))
hist(diff_Glasso)

