library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(RcppProgress)
sourceCpp("./src/ars_multinomial_helper.cpp")
sourceCpp("./src/ars_pois_helper.cpp")

k <- 3
p <- 3
n <- 1000
Z_curr <- matrix(rep(c(1,log(19023),log(19023)),n),n,k,byrow = T)
mu_Z <- matrix(-2.96129,n,k)
Sigma_Z <- 0.357906*diag(k)

y <- matrix(rep(c(0,2300,2300,400),n),nrow = n,byrow = T)

update_Z_helper_Pois(Z_curr,mu_Z,Sigma_Z,y[,-1],k,p,n,100,3,64)
#update_Z_helper_Pois_para(Z_curr,mu_Z,Sigma_Z,y[,-1],k,p,n,700,3,64)
update_Z_helper_multinomial(Z_curr,mu_Z,Sigma_Z,y,k,p,n,100,4,64)
#update_Z_helper_multinomial_para(Z_curr,mu_Z,Sigma_Z,y,k,p,n,1000,10,64)