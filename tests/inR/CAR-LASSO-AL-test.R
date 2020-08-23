library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(ggplot2)
library(reshape2)
library(optimParallel)

rm(list = ls())

get_graph <- function(CAR_sample,k){
    Omega <- matrix(0,k,k)
    Omega[upper.tri(Omega,T)] = apply(CAR_sample$Omega,2,mean)
    Omega <- Omega+t(Omega)
    diag(Omega) <- 0.5 * diag(Omega)
    return(Omega)
}


k <- 5
n_init <- 20 # initial data points
p <- 5

n_step <- 30
n_rep_step <- 5 # number of repeats
n_new <- 4 # number of new experiments each time
n_new_each_step <- n_rep_step * n_new

sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/CAR_active_learning_helper_para.cpp")

# get true parameters
set.seed(12345)
B <- rsparsematrix(k,k,0.2)
omega <- diag(rgamma(k,3,.1))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
Omega <- as.matrix(Omega)

Sigma <- solve(Omega)
beta <- matrix(rnorm(k*p),p,k)
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


res <- data.frame(step = 1:n_step,
    AL_con = NA,
    AL_uncon = NA, rep = NA, random = NA)
X11()
for(i in 1:n_step){

    cat("experiment interation: ",i,"\n")

    # repeat
    cat("  simple repeat:\n")
    design_temp <- design_init[1:n_new_each_step,]
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
    res$rep[i] <- stein_loss_cpp(graph_rep,Omega)

    # random
    cat("  random:\n")
    design_temp <- matrix(rnorm(n_new_each_step*p),ncol = p)
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
    res$rand[i] <- stein_loss_cpp(graph_rand,Omega)


    # AL_con
    cat("  active learning w/ constrain:\n")
    AL_design_con <- optim(par = rnorm(n_new * p),
        fn = expected_G_det_para, 
        old_design = design_AL_con,
        CAR_model = sample_AL_con,
        k = k, p = p, n_new = n_new, n_old = nrow(design_AL_con),
        control = list(maxit = 3000,fnscale = -1),
        method = "L-BFGS-B",lower = -3,upper = 3)
    AL_design_con <- matrix(AL_design_con$par,ncol = p)
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
    res$AL_con[i] <- stein_loss_cpp(graph_AL_con,Omega)


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
    dev.off()
    X11()
    plot(res$step,res$AL_con,type = "l",col = "darkred",ylim = c(0,0.3))
    #lines(res$step,res$AL_uncon)
    lines(res$step,res$rand,col = "blue")
    lines(res$step,res$rep,col = "#058505")
    write.csv(res,"./tests/inR/active_learning_test1.csv")
}
