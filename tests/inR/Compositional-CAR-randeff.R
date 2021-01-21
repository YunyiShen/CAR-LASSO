library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k <- 8
n <- 500
p <- 2
pr <- 2
m <- 2

sourceCpp("./src/Multinomial-CAR-ALASSO-randeff.cpp")

sourceCpp("./src/Multinomial-CAR-ALASSO.cpp")
source("./tests/Formal/Accurancy/Graph_generator.R")

set.seed(12345)
Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
#image(Omega)

Sigma <- Graph_raw$Sigma

Design <- matrix(rnorm(n*p,0,1),n,p)
Design <- (Design-mean(Design))/sd(Design)

beta <- matrix(rnorm(p*k,0,1),p,k)

mu <- rnorm(k,0,1)

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)
Y <- matrix(NA,n,k+1)


design <- matrix(rnorm(n*p),n,p)
design_r <- lapply(1:n,
    function(i,pr){
        w = rep(0,pr);w[sample(pr,1)]=1;return(w)}
    ,pr)

design_r <- Reduce(rbind,design_r)

beta <- matrix(rnorm(p*k),p,k)
mu <- rnorm(k)
xi <- matrix(rgamma(k*m,2,1),m,k)

Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
Sigma <- Graph_raw$Sigma
membership <- matrix(0,pr,m)
membership[1:min(pr,m),1:min(pr,m)] <- diag(min(pr,m))
membership[rowSums(membership)==0,m] <- 1

invsigma_randeff <- membership %*% xi
nu <- matrix(rnorm(k*pr,0,1/sqrt(invsigma_randeff)),pr,k)

Z <- matrix(0,n,k)
Xbeta <- design %*% beta
Xnu <- design_r %*% nu


N <- 500
for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+Xnu[i,]+mu),Sigma)
  p_temp <- c(exp(Z[i,]) , 1)
  Y[i,] <- rmultinom(1, N, p_temp/sum(p_temp))
}

tests_rand <- Multinomial_CAR_ALASSO_randeff_Cpp(data = Y, design = design,
                               design_r = design_r,membership = membership,
                               n_iter = 5000, n_burn_in = 1000, thin_by = 10,
                               r_beta = 1e-2+0*beta, delta_beta = 1e-6 + 0 * beta,
                               r_Omega = rep(1e-2,.5*(k-1)*k),
                               delta_Omega = rep(1e-6,.5*(k-1)*k),
                               alpha = 1e-8, beta = 1e-2,
                               ns = 100,m = 10, emax = 64,
                               progress = T)



test <- Multinomial_CAR_ALASSO_Cpp(Y,  Design, n_iter = 5000, 
                            n_burn_in = 1000, thin_by = 10, 
                            r_beta = 1+0*beta, delta_beta = .01 + 0 * beta,
                            r_Omega = rep(1,.5*(k-1)*k),
                            delta_Omega = rep(.01,.5*(k-1)*k),
                            ns = 100,m = 10, emax = 64,
                            progress = T)

rand_Graph <- 0 * Omega
rand_Graph[upper.tri(rand_Graph,T)] <- apply(tests_rand$Omega,2,mean,na.rm = T)
rand_Graph <- rand_Graph+t(rand_Graph)
diag(rand_Graph) <- 0.5 * diag(rand_Graph)

fix_Graph <- 0 * Omega
fix_Graph[upper.tri(fix_Graph,T)] <- apply(test$Omega,2,mean,na.rm = T)
fix_Graph <- fix_Graph+t(fix_Graph)
diag(fix_Graph) <- 0.5 * diag(fix_Graph)
image(fix_Graph)
