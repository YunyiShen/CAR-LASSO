library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 11
n = 8000
p = 1

sourceCpp("./src/SRG-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")

set.seed(42)
Omega <- rsparsematrix(k,k,0.5)
Omega <- Omega*t(Omega)
diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)

Sigma <- solve(Omega)

Design <- matrix(rnorm(n*p,0,3),n,p)
Design <- (Design-mean(Design))/sd(Design)
#beta <- rsparsematrix(p,k,0.7)
beta <- matrix(rnorm(p*k,1.1,0.05),p,k)
#beta <- as.matrix(beta)
#beta[sample(length(beta),floor(0.3*length(beta)))] <- 0
beta
mu <- rnorm(k)
mu

Xbeta <- Design %*% beta

simu <- matrix(NA,n,k)

for( i in 1:n ){
  simu[i,] <- MASS::mvrnorm(1,Xbeta[i,]+mu,Sigma)
}

test <- SRG_LASSO_Cpp(simu,  Design, n_iter = 2000, 
                      n_burn_in = 1000, thin_by = 10, 
                      r_beta = 1, delta_beta = .1,
                      r_Omega = 1,delta_Omega = .1,
                      progress = T)

Omega_uptri <- Omega[upper.tri(Omega,T)]
diff_SRG <- apply(test$Omega,1,function(w,k){(w-k)/k},Omega_uptri)
hist(diff_SRG)

Glasso <- Graphical_LASSO_Cpp(simu,2000,1000,10,1,.1,T)
diff_Glasso <- apply(Glasso$Omega,1,function(w,k){(w-k)/k},c(Omega))
hist(diff_Glasso)

Glasso_ori <- blockGLasso(simu)

diff_glasso_ori <- sapply(Glasso_ori$Omegas,function(w,Omega){(w[upper.tri(w)]-Omega[upper.tri(w)])/Omega[upper.tri(w)]},Omega)
hist(diff_glasso_ori)
