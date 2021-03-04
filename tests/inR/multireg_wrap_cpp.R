radian

library(Rcpp)
library(RcppProgress)
library(RcppArmadillo)
library(Matrix)

ddd <- CARlasso::simu_AR1(k=5,n=5000)

data1 <- as.matrix(ddd[,1:5])
design <- as.matrix(ddd[,1:5+5])


sourceCpp("./src/multireg_wrap.cpp")
source("./R/multireg-wrap.R")
www <- CAR_multireg_cpp(data1, 
    design, 1000, 
    matrix(0,6,5), 
    diag(1e-8,6,6), 3,
    3*diag(2, 5,5)) 

B_hat <- matrix(colMeans(www$beta), 5,5)
Omega_hat <- CARlasso:::get_graph(www ,5)

www <- CAR_multireg(data1, 
    design, 1000) 




# after Compositional-CAR
mltno <- Multinomial_CAR_multireg_cpp(Y, Design, 10,1000,1,matrix(0,p+1,k),diag(1e-8,p+1,p+1),
3,3 * diag(2,k,k),100,20,64) 
B_hat <- matrix(colMeans(mltno$beta), p,k)
Omega_hat <- CARlasso:::get_graph(mltno ,k)

mltnoR <- Multinomial_CAR_multireg(Y, Design, 10,1000,1,matrix(0,p+1,k),diag(1e-8,p+1,p+1),
3,3 * diag(2,k,k),100,20,64) 

B_hatR <- matrix(colMeans(mltnoR$beta), p,k)
Omega_hatR <- CARlasso:::get_graph(mltnoR ,k)


pois <- Pois_CAR_multireg(Y[,-(k+1)], Design, 10,1000,1) 
B_hatR <- matrix(colMeans(pois$beta), p,k)
Omega_hatR <- CARlasso:::get_graph(pois ,k)


pois <- Pobit_CAR_multireg(1*(Y[,-(k+1)]>1000), Design, 10,1000,1) 

Omega_hatR <- CARlasso:::get_graph(pois ,k)
