library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(dplyr)

rm(list = ls())

sourceCpp("./src/Multinomial-CAR-ALASSO.cpp")
source("./R/multireg-wrap.R")
comp_mat <- read.csv("./real_data/Soil/clean_data/genus_mat_.01_50_without_unclass.csv",row.names = 1) %>% as.matrix()
#comp_mat <- comp_mat[,-36]

Design <- read.csv("real_data/Soil/clean_data/Design_.01_50_without_unclass.csv",row.names = 1)

Design_temp <- Design
Design_temp$res <- 1

Design_dummy <- model.matrix(res~.,data = Design_temp)
Design_dummy <- Design_dummy[,-c(1)]
Design_dummy <- Design_dummy[,-c(2,3,5,7)]
Design_dummy <- apply(Design_dummy,2,function(w){(w-mean(w))/sd(w)})
#Design_dummy <- Design_dummy[,c(2,4,5)]

k <- ncol(comp_mat)-1
p <- ncol(Design_dummy)

Omega <- matrix(0,k,k)

test <- Multinomial_CAR_ALASSO_Cpp(comp_mat,  Design_dummy, n_iter = 50000, 
                                   n_burn_in = 10000, thin_by = 25, 
                                   r_beta = matrix(1e-2,p,k), delta_beta = matrix(1e-6,p,k),
                                   r_Omega = rep(1e-2,.5*(k-1)*k),
                                   delta_Omega = rep(1e-6,.5*(k-1)*k),
                                   ns = 100,m = 10, emax = 64,
                                   progress = T)

A_Graph <- 0 * Omega
A_Graph[upper.tri(A_Graph,T)] <- apply(test$Omega,2,mean,na.rm = T)
A_Graph <- A_Graph+t(A_Graph)
diag(A_Graph) <- 0.5 * diag(A_Graph)
A_beta <- matrix(apply(test$beta,2,mean,na.rm = T),p,k)

multireg_res <- Multinomial_CAR_multireg(data = comp_mat, design = Design_dummy,
                                        n_burn_in = 5000, n_sample = 25000, thin_by = 10)

multireg_Graph <- 0 * Omega
multireg_Graph[upper.tri(multireg_Graph,T)] <- apply(multireg_res$Omega,2,mean,na.rm = T)
multireg_Graph <- multireg_Graph+t(multireg_Graph)
diag(multireg_Graph) <- 0.5 * diag(multireg_Graph)
multireg_beta <- matrix(apply(multireg_res$beta,2,mean,na.rm = T),p,k)

Graph_binary <- abs(A_Graph/multireg_Graph) > .5
beta_binary <- abs(A_beta/multireg_beta) > .5


save.image("./real_data/Soil/res/CAR_full_design_genus_.01_50_without_unclass_lambda_prior_1e-2_1e-6_long_chain.RData")  


