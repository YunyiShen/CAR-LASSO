library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(ggplot2)
library(reshape2)
library(optimParallel)
library(RcppParallel)

rm(list = ls())
setThreadOptions(numThreads = defaultNumThreads() - 2)
get_graph <- function(CAR_sample,k){
  Omega <- matrix(0,k,k)
  Omega[upper.tri(Omega,T)] = apply(CAR_sample$Omega,2,mean)
  Omega <- Omega+t(Omega)
  diag(Omega) <- 0.5 * diag(Omega)
  return(Omega)
}

log_l2_upt_dist <- function(Omega_hat,Omega){
  log(sum((Omega_hat[upper.tri(Omega,T)] -
         Omega[upper.tri(Omega,T)])^2))
  
}

G_det_rep <- function(new_design,
                      old_design,
                      CAR_model,
                      k, p, n_new, n_old,
                      nrep,pan_scale=.1, pan_amp=.01){
  temp <- matrix(new_design,ncol = p)
  temp <- lapply(1:nrep,function(i){temp})
  temp <- Reduce(rbind,temp)
  temp <- c(temp)
  expected_G_det_para(temp,old_design,
                      CAR_model,
                      k, p, n_new * nrep, n_old) - 
    sum(pan_amp * exp(pan_scale * new_design^2))
  
  
}

G_tr_rep <- function(new_design,
                      old_design,
                      CAR_model,
                      k, p, n_new, n_old,
                      nrep,pan_scale=.1, pan_amp=.01){
  temp <- matrix(new_design,ncol = p)
  temp <- lapply(1:nrep,function(i){temp})
  temp <- Reduce(rbind,temp)
  temp <- c(temp)
  expected_G_tr_para(temp,old_design,
                      CAR_model,
                      k, p, n_new * nrep, n_old) - 
    sum(pan_amp * exp(pan_scale * new_design^2))
  
  
}



full_tr_rep <- function(new_design,
                     old_design,
                     CAR_model,
                     k, p, n_new, n_old,
                     nrep,pan_scale=.1, pan_amp=.01){
  temp <- matrix(new_design,ncol = p)
  temp <- lapply(1:nrep,function(i){temp})
  temp <- Reduce(rbind,temp)
  temp <- c(temp)
  expected_full_tr_para(temp,old_design,
                     CAR_model,
                     k, p, n_new * nrep, n_old) - 
    sum(pan_amp * exp(pan_scale * new_design^2))
  
  
}

k <- 10
n_init <- 30 # initial data points
p <- 5

n_step <- 30
n_rep_step <- 2 # number of repeats
n_new <- 5 # number of new experiments each time
n_new_each_step <- n_rep_step * n_new

sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/CAR_active_learning_helper_para.cpp")

# get true parameters
set.seed(12345)
B <- rsparsematrix(k,k,0.2)
omega <- diag(rgamma(k,2,.1))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
Omega <- as.matrix(Omega)

Sigma <- solve(Omega)
beta <- matrix(rnorm(k*p,0,2),p,k)
mu <- matrix(rnorm(k))


design_init <- matrix(rnorm(p * n_init), ncol = p)
Xbeta_init <- design_init %*% beta

design_AL_con <- design_init
design_AL_uncon <- design_init
design_rand <- design_init
design_rep <- design_init

data_init <- matrix(NA,n_init,k)
for( i in 1:n_init ){
  data_init[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta_init[i,]+mu),Sigma)
}

data_AL_con <- data_init
data_AL_uncon <- data_init
data_rand <- data_init
data_rep <- data_init

data_temp <- matrix(NA,n_new_each_step,k)


# initial result
sample_init <- CAR_LASSO_Cpp(data_init,  design_init, n_iter = 25000, 
                             n_burn_in = 5000, thin_by = 50, 
                             r_beta = 1, delta_beta = .01,
                             r_Omega = 1,delta_Omega = .01,
                             progress = T)
sample_AL_con <- sample_AL_uncon <- sample_rand <- sample_rep <- sample_init

graph_init <- get_graph(sample_init,k)

res1 <- data.frame(step = 0:n_step,
                  AL_con = NA,
                  AL_uncon = NA, rep = NA, rand = NA)
res2 <- data.frame(step = 0:n_step,
                   AL_con = NA,
                   AL_uncon = NA, rep = NA, rand = NA)
res1[1,2:5] <- stein_loss_cpp(graph_init,Omega)
res2[1,2:5] <- log_l2_upt_dist(graph_init,Omega)
par(mfrow = c(1,2))
#X11()
for(i in 1:n_step+1){
  
  cat("experiment interation: ",i-1,"\n")
  
  # repeat
  cat("  simple repeat:\n")
  design_temp <- design_init[1,]
  design_temp <- lapply(1:n_new_each_step,function(i){design_temp})
  design_temp <- Reduce(rbind, design_temp)
  Xbeta_temp <- design_temp %*% beta
  for( j in 1:n_new_each_step ){
    data_temp[j,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta_temp[j,]+mu),Sigma)
  }
  data_rep <- rbind(data_rep,data_temp)
  design_rep <- rbind(design_rep,design_temp)
  
  sample_rep <- CAR_LASSO_Cpp(data_rep,  design_rep, n_iter = 25000, 
                              n_burn_in = 5000, thin_by = 10, 
                              r_beta = 1, delta_beta = .01,
                              r_Omega = 1,delta_Omega = .01,
                              progress = T)
  
  graph_rep <- get_graph(sample_rep,k)
  res1$rep[i] <- stein_loss_cpp(graph_rep,Omega)
  res2$rep[i] <- log_l2_upt_dist(graph_rep,Omega)
  
  # random
  cat("  random:\n")
  design_temp <- matrix(rnorm(n_new*p),ncol = p)
  design_temp <- lapply(1:n_rep_step,function(i){design_temp})
  design_temp <- Reduce(rbind, design_temp)
  Xbeta_temp <- design_temp %*% beta
  for( j in 1:n_new_each_step ){
    data_temp[j,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta_temp[j,]+mu),Sigma)
  }
  data_rand <- rbind(data_rand,data_temp)
  design_rand <- rbind(design_rand,design_temp)
  
  sample_rand <- CAR_LASSO_Cpp(data_rand,  design_rand, n_iter = 25000, 
                               n_burn_in = 5000, thin_by = 50, 
                               r_beta = 1, delta_beta = .01,
                               r_Omega = 1,delta_Omega = .01,
                               progress = T)
  
  graph_rand <- get_graph(sample_rand,k)
  res1$rand[i] <- stein_loss_cpp(graph_rand,Omega)
  res2$rand[i] <- log_l2_upt_dist(graph_rand,Omega)
  
  
  # AL_con
  cat("  active learning w/ constrain:\n")
  AL_design_con <- optim(par = rep(0,n_new * p),
                         fn = G_tr_rep, 
                         old_design = design_AL_con,
                         CAR_model = sample_AL_con,
                         k = k, p = p, n_new = n_new, n_old = nrow(design_AL_con),
                         nrep = n_rep_step,
                         control = list(maxit = 3000,fnscale = -1)
                         ,method = "L-BFGS-B"
                         #,lower = -3,upper = 3
                         #,method = "SANN"
                         )
  AL_design_con <- matrix(AL_design_con$par,ncol = p)
  cat(AL_design_con,"\n")
  design_temp <- lapply(1:n_rep_step,function(i){AL_design_con})
  design_temp <- Reduce(rbind, design_temp)
  Xbeta_temp <- design_temp %*% beta
  
  for( j in 1:n_new_each_step ){
    data_temp[j,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta_temp[j,]+mu),Sigma)
  }
  data_AL_con <- rbind(data_AL_con,data_temp)
  design_AL_con <- rbind(design_AL_con,design_temp)
  
  sample_AL_con <- CAR_LASSO_Cpp(data_AL_con,  design_AL_con, n_iter = 25000, 
                                 n_burn_in = 5000, thin_by = 50, 
                                 r_beta = 1, delta_beta = .01,
                                 r_Omega = 1,delta_Omega = .01,
                                 progress = T)
  
  graph_AL_con <- get_graph(sample_AL_con,k)
  res1$AL_con[i] <- stein_loss_cpp(graph_AL_con,Omega)
  res2$AL_con[i] <- log_l2_upt_dist(graph_AL_con,Omega)
  
  
  # AL_uncon
  #cat("  active learning w/o constrain:\n")
  #AL_design_uncon <- optim(par = rnorm(n_new * p),
  #    fn = expected_G_det, 
  #    old_design = design_AL_uncon,
  #    CAR_model = sample_AL_uncon,
  #    k = k, p = p, n_new = n_new, n_old = nrow(design_AL_uncon),
  #    control = list(maxit = 3000,fnscale = -1),
  #    method = "L-BFGS-B")
  #AL_design_uncon <- matrix(AL_design_uncon$par,ncol = p)
  #design_temp <- lapply(1:n_rep_step,function(i){AL_design_uncon})
  #design_temp <- Reduce(rbind, design_temp)
  #Xbeta_temp <- design_temp %*% beta
  
  #for( j in 1:n_new_each_step ){
  #    data_temp[j,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta_temp[j,]+mu),Sigma)
  #}
  #data_AL_uncon <- rbind(data_AL_uncon,data_temp)
  #design_AL_uncon <- rbind(design_AL_uncon,design_temp)
  
  #sample_AL_uncon <- CAR_LASSO_Cpp(data_AL_uncon,  design_AL_uncon, n_iter = 25000, 
  #                      n_burn_in = 5000, thin_by = 50, 
  #                      r_beta = 1, delta_beta = .01,
  #                      r_Omega = 1,delta_Omega = .01,
  #                      progress = T)
  
  #graph_AL_uncon <- get_graph(sample_AL_uncon,k)
  #res$AL_uncon[i] <- stein_loss_cpp(graph_AL_uncon,Omega)
  #dev.off()
  #X11()
  plot(res1$step,res1$AL_con,type = "l",col = "darkred",ylim = c(0,max(res1[,2:5],na.rm = T)))
  #lines(res$step,res$AL_uncon)
  lines(res1$step,res1$rand,col = "blue")
  lines(res1$step,res1$rep,col = "#058505")
  
  plot(res2$step,res2$AL_con,type = "l",col = "darkred",ylim = c(0,max(res2[,2:5],na.rm = T)))
  #lines(res$step,res$AL_uncon)
  lines(res2$step,res2$rand,col = "blue")
  lines(res2$step,res2$rep,col = "#058505")
  write.csv(res1,"./tests/inR/active_learning_test1.csv")
  write.csv(res2,"./tests/inR/active_learning_test1.csv")
}
