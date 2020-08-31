library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(ggplot2)
library(reshape2)
library(optimParallel)
library(RcppParallel)
library(BB)
library(GenSA)

sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/CAR_active_learning_helper_para.cpp")
sourceCpp("./src/CAR_FI_helper.cpp")

get_graph <- function(CAR_sample,k){
  Omega <- matrix(0,k,k)
  Omega[upper.tri(Omega,T)] = apply(CAR_sample$Omega,2,mean)
  Omega <- Omega+t(Omega)
  diag(Omega) <- 0.5 * diag(Omega)
  return(Omega)
}


Prior_Info_det <- function(design_vec,CAR_sample,k,pp,nrep,pan_scale=.1,pan_amp=.01){
    p = pp
    Omega_hat <- get_graph(CAR_sample,k)
    beta_hat <- matrix(colMeans(CAR_sample$beta),nrow = p)
    mu_hat <- colMeans(CAR_sample$mu)
    Emp_cov_prior <- cov(Reduce(cbind,CAR_sample[1:3]))

    Emp_info_prior <- solve(Emp_cov_prior)

    designmat <- matrix(design_vec,ncol = p)
    designmat <- lapply(1:nrep,function(i) designmat)
    designmat <- Reduce(rbind,designmat)

    infomat <- CAR_FI_para(designmat,Omega_hat,beta_hat,mu_hat,k,p=p) + Emp_info_prior
    return(-log(det(infomat[1:((p+1)*k),1:((p+1)*k)]))+ 
        sum(pan_amp * exp(pan_scale * designmat^2))/length(designmat))
}

Prior_Info_tr <- function(design_vec,CAR_sample,k,pp,nrep,pan_scale=.1,pan_amp=.01){
  p = pp
  Omega_hat <- get_graph(CAR_sample,k)
  beta_hat <- matrix(colMeans(CAR_sample$beta),nrow = p)
  mu_hat <- colMeans(CAR_sample$mu)
  Emp_cov_prior <- cov(Reduce(cbind,CAR_sample[1:3]))
  
  Emp_info_prior <- solve(Emp_cov_prior)
  
  designmat <- matrix(design_vec,ncol = p)
  designmat <- lapply(1:nrep,function(i) designmat)
  designmat <- Reduce(rbind,designmat)
  
  infomat <- CAR_FI_para(designmat,Omega_hat,beta_hat,mu_hat,k,p=p) + Emp_info_prior
  return(log(sum(diag(solve(infomat))))+ 
           sum(pan_amp * exp(pan_scale * designmat^2))/length(designmat))
}



log_l2_upt_dist <- function(Omega_hat,Omega){
  log(sum((Omega_hat[upper.tri(Omega)] -
         Omega[upper.tri(Omega)])^2))
  
}

log_l2_dist <- function(beta_hat,beta){
  log(sum(
    (beta_hat-beta)^2
  ))
  
}


k <- 12
n_init <- 200 # initial data points
p <- 12

n_step <- 15
n_rep_step <- 1 # number of repeats
n_new <- 30 # number of new experiments each time
n_new_each_step <- n_rep_step * n_new


# get true parameters
#set.seed(12345)
deg = T
while(deg){
  B <- rsparsematrix(k,k,0.2)
  omega <- diag(rgamma(k,5,.1))
  I <- diag(rep(1,k))
  Omega <- t(I-B) %*% omega %*% (I-B)
  Omega <- as.matrix(Omega)
  Sigma <- solve(Omega)
  deg = kappa(Sigma)>1e4
}
cat(kappa(Sigma),"\n")


#beta <- matrix(rnorm(k*p,0,2),p,k)
beta <- as.matrix(rsparsematrix(p,k,0.5))
mu <- matrix(rnorm(k,0,1))


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
#res1[1,2:5] <- stein_loss_cpp(graph_init,Omega)
res1[1,2:5] <- log_l2_dist(colMeans(sample_init$beta),c(beta))
#res2[1,2:5] <- log_l2_dist(graph_init[upper.tri(Omega,T)],
#                           Omega[upper.tri(Omega,T)])
res2[1,2:5] <- stein_loss_cpp(graph_init,Omega)
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
                              n_burn_in = 5000, thin_by = 50, 
                              r_beta = 1, delta_beta = .01,
                              r_Omega = 1,delta_Omega = .01,
                              progress = T)
  
  graph_rep <- get_graph(sample_rep,k)
  res1$rep[i] <- log_l2_dist(colMeans(sample_rep$beta),c(beta))
    #stein_loss_cpp(graph_rep,Omega)
  #res2$rep[i] <- log_l2_dist(graph_rep[upper.tri(Omega,T)],
  #           Omega[upper.tri(Omega,T)])
  res2$rep[i] <- stein_loss_cpp(graph_rep,Omega)
  
  # random
  cat("  random:\n")
  design_temp <- matrix(rnorm(n_new*p,0,1),ncol = p)
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
  res1$rand[i] <- log_l2_dist(colMeans(sample_rand$beta),c(beta))
    #stein_loss_cpp(graph_rand,Omega)
  #res2$rand[i] <- log_l2_dist(graph_rand[upper.tri(Omega,T)],
  #                           Omega[upper.tri(Omega,T)])
  res2$rand[i] <- stein_loss_cpp(graph_rand,Omega)
  
  
  # AL_con
  cat("  active learning:\n")
  guess <- rnorm(n_new * p,0,1)
  AL_design_con <- BBoptim(par = guess,
                         fn = Prior_Info_tr,
                         #lower = -3 + 0*guess,upper = 3+0*guess,
                         CAR_sample = sample_AL_con,
                         k = k, 
                         pp = p,
                         nrep = n_rep_step,
                         pan_amp = 1*5e-3,
                        #control = list(max.call = 1e5,
                         #               verbose = T,
                          #             nb.stop.improvement = 1e-5,
                           #             simple.function = T)
                         #,method = "L-BFGS-B"
                         #,method = 1
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
  res1$AL_con[i] <- log_l2_dist(colMeans(sample_AL_con$beta),c(beta))
    #stein_loss_cpp(graph_AL_con,Omega)
  #res2$AL_con[i] <- log_l2_dist(graph_AL_con[upper.tri(Omega,T)],
  #                           Omega[upper.tri(Omega,T)])
  res2$AL_con[i] <- stein_loss_cpp(graph_AL_con,Omega)
  
  
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
  plot(res1$step,res1$AL_con,type = "l",col = "darkred",
       ylim = c(min(res1[,2:5]-.25,na.rm = T),
                max(res1[,2:5]+.25,na.rm = T)),
       xlab = "steps",
       ylab = "log L2 loss beta")
  #lines(res$step,res$AL_uncon)
  lines(res1$step,res1$rand,col = "blue")
  lines(res1$step,res1$rep,col = "#058505")
  
  #legend("topright", 
  #       legend=c("AL", "Rand" , "Rep"),
  #       col=c("darkred", "blue" , "#058505"), 
  #       lty=1, box.lty=0)
  
  plot(res2$step,res2$AL_con,type = "l",col = "darkred",
       ylim = c(0,max(res2[,2:5],na.rm = T)),
       xlab = "steps",
       ylab = "stein's loss of Omega")
  #lines(res$step,res$AL_uncon)
  lines(res2$step,res2$rand,col = "blue")
  lines(res2$step,res2$rep,col = "#058505")
  #legend("topright", cex=.5,
  #       legend=c("AL", "Rand" , "Rep"),
  #       col=c("darkred", "blue" , "#058505"), 
  #       lty=1, box.lty=1)
  write.csv(res1,"./tests/inR/active_learning_test1.csv")
  write.csv(res2,"./tests/inR/active_learning_test2.csv")
}
