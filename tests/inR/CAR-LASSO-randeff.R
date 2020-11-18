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


sourceCpp("./src/CAR-LASSO-randeff.cpp")



tests <- CAR_LASSO_randeff_Cpp( data = data1,     
                            design = design,   
                            design_r = design_r,
                            membership = membership, 
                           n_iter = 10,          
                           n_burn_in = 10,       
                           thin_by = 1,         
                            r_beta = 1,       
                            delta_beta = 0.01,
                            r_Omega = 1, 
                            delta_Omega = 0.01,
                            alpha = 0.01, 
                            beta = 0.01,  
                           progress = T)
