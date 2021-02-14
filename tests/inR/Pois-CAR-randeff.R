library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k <- 5
n <- 500
p <- 8
pr <- 10
m <- 1

sourceCpp("./src/Pois-CAR-ALASSO-randeff.cpp")
sourceCpp("./src/Pois-CAR-ALASSO.cpp")
sourceCpp("./src/Pois-CAR-LASSO.cpp")
source("./tests/Formal/Accurancy/Graph_generator.R")

set.seed(12345)

design <- matrix(rnorm(n*p),n,p)
design_r <- lapply(1:n,
    function(i,pr){
        w = rep(0,pr);w[sample(pr,1)]=1;return(w)}
    ,pr)

design_r <- Reduce(rbind,design_r)

beta <- matrix(rnorm(p*k,0,.1),p,k)
mu <- rnorm(k,0,1)
xi <- matrix(rgamma(k*m,10,1),m,k)

Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
Sigma <- Graph_raw$Sigma

#Sigma <- 0.1 * diag(k)
membership <- matrix(0,pr,m)
membership[1:min(pr,m),1:min(pr,m)] <- diag(min(pr,m))
membership[rowSums(membership)==0,m] <- 1

invsigma_randeff <- membership %*% xi
nu <- matrix(rnorm(k*pr,0,1/sqrt(invsigma_randeff)),pr,k)

Y <- Z <- matrix(0,n,k)
Xbeta <- design %*% beta
Xnu <- design_r %*% nu



for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+Xnu[i,]+mu),Sigma)
  Y[i,] <- rpois(k,exp(Z[i,]))
}

hist(exp(Z))
hist(Y)

tests_rand <- Pois_CAR_ALASSO_randeff_Cpp(data = Y, design = design,
                               design_r = design_r,membership = membership,
                               n_iter = 5000, n_burn_in = 1000, thin_by = 10,
                               r_beta = 1+0*beta, delta_beta = 1e-2 + 0 * beta,
                               r_Omega = rep(1,.5*(k-1)*k),
                               delta_Omega = rep(1e-2,.5*(k-1)*k),
                               alpha = 0.01, beta = 0.01,
                               ns = 100,m = 10, emax = 64,
                               progress = T)

rand_Graph <- 0 * Omega
rand_Graph[upper.tri(rand_Graph,T)] <- apply(tests_rand$Omega,2,mean,na.rm = T)
rand_Graph <- rand_Graph+t(rand_Graph)
diag(rand_Graph) <- 0.5 * diag(rand_Graph)

tests_fix <- Pois_CAR_ALASSO_Cpp(data = Y, design = design,
                               #design_r = design_r,membership = membership,
                               n_iter = 5000, n_burn_in = 1000, thin_by = 10,
                               r_beta = 1e-2+0*beta, delta_beta = 1e-6 + 0 * beta,
                               r_Omega = rep(1e-2,.5*(k-1)*k),
                               delta_Omega = rep(1e-6,.5*(k-1)*k),
                               #alpha = 1, beta = 1,
                               ns = 100,m = 10, emax = 64,
                               progress = T)

tests_fix <- Pois_CAR_LASSO_Cpp(data = Y, design = design,
                               #design_r = design_r,membership = membership,
                               n_iter = 5000, n_burn_in = 1000, thin_by = 10,
                               r_beta = 1, delta_beta = 1e-2,
                               r_Omega = 1,
                               delta_Omega = 1e-2,
                               #alpha = 1, beta = 1,
                               ns = 100,m = 10, emax = 64,
                               progress = T)
