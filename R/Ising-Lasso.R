library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
require(coda)

rm(list = ls())
sourceCpp("./src/Ising-LASSO.cpp")

n = 800 # number of data
p = 1 # number of predictors
k = 5 # number of nodes

set.seed(42)
graph <- matrix(rnorm(k^2,0,0.1),k,k)
diag(graph) <- 0
graph <- graph+t(graph)

zrs <- sample(1:k,4)
zrs_ind <- sample(4,2)
zrs_x <- c(zrs[zrs_ind],zrs[-zrs_ind])
zrs_y <- c(zrs[-zrs_ind],zrs[zrs_ind])
graph[cbind(zrs_x,zrs_y)] <- 0
isSymmetric(graph)

beta <-  matrix( rnorm(k*p),p,k)
mu <- 0.5*rnorm(k)
design <- matrix(rnorm(n*p),n,p)

thr_mat <- design %*% beta 
thr_mat <- apply(thr_mat,1,`+`,mu)
#thr_mat <- t(thr_mat)

test_data <- matrix(NA,k,n)

for(i in 1:n){
  test_data[,i] <- IsingSamplerCpp(1,graph,thr_mat[,i],100,c(-1,1),T)
}
cor_emp <- cor(t(test_data))
par(mfrow = c(1,2))
image(graph[,k-1:k+1])
image(cor_emp[,k-1:k+1])

test <- Ising_LASSO_Cpp(data = test_data, design = design, n_iter = 15000, 
                        n_burn_in = 5000, thin_by = 10, r_beta = 2, 
                        delta_beta = .5, r_J = 2, delta_J = .5, 
                        propsd_mu = rep(0.025,k), propsd_beta = rep(0.025,k*p), 
                        propsd_J = rep(0.01,0.5*(k*(k-1))), 
                        propsd_lambda = c(.5,.5), 
                        exact = TRUE, nauxIter=100, progress=F, verbos=T, 
                        reportby=100)


