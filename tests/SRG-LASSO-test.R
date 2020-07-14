library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 11
n = 1100
p = 1

sourceCpp("./src/Probit-SRG-LASSO.cpp")
sourceCpp("./src/Probit-Graphical-LASSO.cpp")
sourceCpp("./src/SRG-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")
source("./R/SRG-LASSO.R")


#set.seed(42)
B <- rsparsematrix(k,k,0.3)
omega <- diag(rgamma(k,3,.1))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
#diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)
image(Omega)

Sigma <- solve(Omega)

Design <- matrix(rnorm(n*p,0,1),n,p)
Design <- (Design-mean(Design))/sd(Design)
colnames(Design) <- paste0("x",1:p)


beta <- matrix(rnorm(p*k,2,1),p,k)
#beta[sample(p*k,floor(0.3*p*k))] = 0

mu <- rnorm(k)
#mu

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)
Y <- Z

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Xbeta[i,]+mu,Sigma)
  Y[i,] <- 1 * ((rnorm(k,Z[i,],1))>0)
}
colnames(Y) <- paste0("y",1:k)
colnames(Z) <- paste0("z",1:k)

full_data <- data.frame(Y,Z,Design)

formula_Z <- paste0(paste0("z",1:k,collapse = "+"),"~",paste0("x",1:p,collapse = "+"))
formula_Z <- as.formula(formula_Z)

test_R <- SRGlasso(formula_Z,full_data,verbos = T,n_iter = 8000,n_burn_in = 1000)


test <- Proit_SRG_LASSO_Cpp(Y,  Design, n_iter = 2000, 
                            n_burn_in = 1000, thin_by = 10, 
                            r_beta = 1, delta_beta = .1,
                            r_Omega = 1,delta_Omega = .1,
                            progress = T)

par(mfrow = c(1,2))
Omega_uptri <- Omega[upper.tri(Omega,T)]
diff_probit <- apply(test$Omega,1,function(w,k){(w-k)/w},Omega_uptri)
hist(diff_probit)

SRG_test <- SRG_LASSO_Cpp(Z,  Design, n_iter = 5000, 
                          n_burn_in = 1000, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)

diff_SRG <- apply(SRG_test$Omega,1,function(w,k){(w-k)},Omega_uptri)
hist(diff_SRG)

SRG_Graph <- 0 * Omega
SRG_Graph[upper.tri(SRG_Graph,T)] = apply(SRG_test$Omega,2,median)
SRG_Graph = SRG_Graph+t(SRG_Graph)
diag(SRG_Graph) = 0.5 * diag(SRG_Graph)
image((SRG_Graph)>1*sd(SRG_Graph[upper.tri(Omega)]))
hist((SRG_Graph-Omega)/Omega)
image(Omega>1*sd(Omega[upper.tri(Omega)]))

probit_Glasso <-Probit_Graphical_LASSO_Cpp(Y,2000,1000,10,1,.1,T)
diff_pGlasso <- apply(probit_Glasso$Omega,1,function(w,k){(w-k)},Omega_uptri)

hist(diff_pGlasso)




Glasso <- Graphical_LASSO_Cpp(Z,5000,1000,10,1,.01,T)
diff_Glasso <- apply(Glasso$Omega,1,function(w,k){(w-k)/k},Omega_uptri)

hist(diff_Glasso)

Glasso_Graph <- 0 * Omega
Glasso_Graph[upper.tri(Glasso_Graph,T)] = apply(Glasso$Omega,2,median)
Glasso_Graph = Glasso_Graph+t(Glasso_Graph)
diag(Glasso_Graph) = 0.5 * diag(Glasso_Graph)
image((Glasso_Graph>sd(Glasso_Graph[upper.tri(Omega)])))
hist((Glasso_Graph-Omega)/Omega)
image(Omega)

