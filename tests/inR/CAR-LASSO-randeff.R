library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k <- 3
n <- 30
p <- 2
pr <- 2
m <- 2

sourceCpp("./src/CAR_LASSO_randeff_helper.cpp")

data1 <- matrix(rnorm(n*k),n,k)
centered_data <- 0 * data1
design <- matrix(rnorm(n*p),n,p)
design_r <- (matrix(rnorm(n*pr),n,pr) > 0) + 0.0
beta <- matrix(rnorm(p*k),p,k)
mu <- rnorm(k)
xi <- matrix(rgamma(k*m,1,1),m,k)
Omega <- diag(k)
membership <- matrix(0,pr,m)
membership[1:min(pr,m),1:min(pr,m)] <- diag(min(pr,m))
membership[rowSums(membership)==0,m] <- 1

xi_curr <- matrix( rgamma(k*m,1,1),m,k)

nu <- update_car_nu_helper(data1, design,  
    design_r, membership, beta,mu, xi_curr, Omega, k,  pr,  n)

update_xi_helper(xi_curr, nu, membership,1.0,1.0, k, pr, 1)


get_data_centered(centered_data,data1, design_r, nu,Omega)

## start the overall test for the sampling algorithm
sourceCpp("./src/CAR-LASSO-randeff.cpp")
sourceCpp("./src/CAR-ALASSO-randeff.cpp")
sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/CAR-ALASSO.cpp")
source("./tests/Formal/Accurancy/Graph_generator.R")

set.seed(12345)
k <- 10
n <- 500
p <- 2
pr <- 5
m <- 1

design <- matrix(rnorm(n*p),n,p)
design_r <- lapply(1:n,
    function(i,pr){
        w = rep(0,pr);w[sample(pr,1)]=1;return(w)}
    ,pr)

design_r <- Reduce(rbind,design_r)

beta <- matrix(rnorm(p*k),p,k)
mu <- rnorm(k)
xi <- matrix(rgamma(k*m,1,1),m,k)

Graph_raw <- g_model1(k)
Omega <- Graph_raw$Omega
Sigma <- Graph_raw$Sigma
membership <- matrix(0,pr,m)
membership[1:min(pr,m),1:min(pr,m)] <- diag(min(pr,m))
membership[rowSums(membership)==0,m] <- 1

invsigma_randeff <- membership %*% xi
nu <- matrix(rnorm(k*pr,0,1/sqrt(invsigma_randeff)),pr,k)

Y <- matrix(0,n,k)
Xbeta <- design %*% beta
Xnu <- design_r %*% nu

## generate data
for(i in 1:n){
    Y[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+Xnu[i,]+mu),Sigma)
}


tests_rand <- CAR_LASSO_randeff_Cpp(data = Y, design = design,
                               design_r = design_r,membership = membership,
                               n_iter = 10000, n_burn_in = 1000, thin_by = 10,
                               r_beta = 1, delta_beta = 0.01,
                               r_Omega = 1, delta_Omega = 0.01,
                               alpha = 1e-8, beta = 1e-2,
                               progress = T)


CAR_rand_Graph <- 0 * Omega
CAR_rand_Graph[upper.tri(CAR_rand_Graph,T)] <- apply(tests_rand$Omega,2,mean)
CAR_rand_Graph <- CAR_rand_Graph+t(CAR_rand_Graph)
diag(CAR_rand_Graph) <- 0.5 * diag(CAR_rand_Graph)
image((CAR_rand_Graph))


tests_fix <- CAR_LASSO_Cpp(data = Y, design = design,
                               n_iter = 10000, n_burn_in = 1000, thin_by = 10,
                               r_beta = 1, delta_beta = 0.01,
                               r_Omega = 1, delta_Omega = 0.01,
                               progress = T)


CAR_fix_Graph <- 0 * Omega
CAR_fix_Graph[upper.tri(CAR_fix_Graph,T)] <- apply(tests_fix$Omega,2,mean)
CAR_fix_Graph <- CAR_fix_Graph+t(CAR_fix_Graph)
diag(CAR_fix_Graph) <- 0.5 * diag(CAR_fix_Graph)
image((CAR_fix_Graph))


tests_fix_A <- CAR_ALASSO_Cpp(data = Y, design = design,
                               n_iter = 10000, n_burn_in = 1000, thin_by = 10,
                               r_beta = 1e-2+0*beta, delta_beta = 1e-6 + 0 * beta,
                               r_Omega = rep(1e-2,.5*(k-1)*k),
                               delta_Omega = rep(1e-6,.5*(k-1)*k),
                               lambda_diag = rep(0,k),
                               progress = T)

CAR_fix_A_Graph <- CARlasso:::get_graph(tests_fix_A,k)



stein_loss_cpp(CAR_fix_Graph,Omega)
stein_loss_cpp(CAR_rand_Graph,Omega)
stein_loss_cpp(CAR_fix_A_Graph,Omega)
stein_loss_cpp(CAR_rand_A_Graph,Omega)


# multireg

sourceCpp("./src/multireg-randeff.cpp")

test_multireg <- CAR_multireg_randeff_cpp(data = Y, design = design,
                               design_r = design_r,membership = membership,
                               n_iter = 10000, n_burn_in = 1000, thin_by = 10,
                               Bbar = matrix(0,p+1,k), A = diag(1e-8,p+1,p+1),
                               nu = 3, V = 3*diag(2,k,k), 
                               alpha = 1e-8, beta = 1e-2)

Graph_multireg <- CARlasso:::get_graph(test_multireg, k=k)
CARlasso:::stein_loss_cpp(Graph_multireg, Omega)



tests_randA <- CAR_ALASSO_randeff_Cpp(data = Y, design = design,
                               design_r = design_r,membership = membership,
                               n_iter = 10000, n_burn_in = 1000, thin_by = 10,
                               r_beta = 1e-2+0*beta, delta_beta = 1e-6 + 0 * beta,
                               r_Omega = rep(1e-2,.5*(k-1)*k),
                               delta_Omega = rep(1e-6,.5*(k-1)*k),
                               alpha = 1e-8, beta = 1e-2,
                               lambda_diag = matrix(0, k , 1),
                               progress = T)

CAR_rand_A_Graph <- CARlasso:::get_graph(tests_randA, k=k)
CARlasso:::stein_loss_cpp(CAR_rand_A_Graph, Omega)
