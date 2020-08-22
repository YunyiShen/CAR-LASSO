library(Rcpp)
library(RcppArmadillo)
library(Matrix)
sourceCpp("./src/CAR_FI_helper.cpp")



min_eig_graph_batch <- function(design_vec,FI_pre,Omega,beta,mu,k,p){
    design <- matrix(design_vec,ncol = p)
    design_list <- lapply(1:nrow(design),function(i,w){t(w[i,])},design)

    temp <- lapply(design_list,CAR_FI,Omega,beta,mu,k,p)
    temp <- Reduce("+",temp)
    temp <- temp + FI_pre
    return(-min(abs(eigen(
            temp[(k ) * (p+1) + 1:(.5 * (k + 1) * k),
                 (k ) * (p+1) + 1:(.5 * (k + 1) * k)]
        )$value)))
}

min_eig_all_batch <- function(design_vec,FI_pre,Omega,beta,mu,k,p){
    design <- matrix(design_vec,ncol = p)
    design_list <- lapply(1:nrow(design),function(i,w){t(w[i,])},design)

    temp <- lapply(design_list,CAR_FI,Omega,beta,mu,k,p)
    temp <- Reduce("+",temp)
    temp <- temp + FI_pre
    return(-min(abs(eigen(
            temp
        )$value)))
}

FI_mat_batch <- function(design_vec,FI_pre,Omega,beta,mu,k,p){
    design <- matrix(design_vec,ncol = p)
    design_list <- lapply(1:nrow(design),function(i,w){t(w[i,])},design)

    temp <- lapply(design_list,CAR_FI,Omega,beta,mu,k,p)
    temp <- Reduce("+",temp)
    temp <- temp + FI_pre
    return(temp)
}


k <- 5
p <- 5
batch <- 10
max_exp <- 200

beta <- matrix(rnorm(k * p,0,1) , p,k)
mu <- rnorm(k)

B <- rsparsematrix(k,k,0.1)
omega <- diag(rgamma(k,1,2))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
#Omega <- omega
#diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)

mineig_test <- sapply(1:1000,function(dummy){min_eig_graph_batch(rnorm(5*p),0 * pre_info_AL,Omega,beta,mu,k,p)}) 
hist(mineig_test)


init_batch <- batch
set.seed(42)
init_design <- matrix(rnorm(p*init_batch),init_batch,p)
init_info <- CAR_FI(t(init_design[1,]),Omega,beta ,mu,k,p)
for(i in 2:init_batch){
    init_info <- CAR_FI(t(init_design[i,]),Omega,beta ,mu,k,p) + init_info
}
pre_info_AL <- init_info
pre_info_rep <- init_info
pre_info_rand <- init_info

AL_round <- max_exp/batch
res <- data.frame(n = 1:AL_round,rep = NA, random = NA, AL = NA)

for(i in 1:AL_round){
    pre_info_rep <- FI_mat_batch(c(init_design), 
                                  pre_info_rep,
                                  Omega,beta,mu,k,p)
    res$rep[i] = min(abs(eigen(
        pre_info_rep[(k ) * (p+1) + 1:(.5 * (k + 1) * k),
        k * (p+1) + 1:(.5 * (k + 1) * k)])$value))

    pre_info_rand <- FI_mat_batch(rnorm(batch * p), 
                                  pre_info_rand,
                                  Omega,beta,mu,k,p)
    res$random[i] = min(abs(eigen(
        pre_info_rand[(k ) * (p+1) + 1:(.5 * (k + 1) * k),
        k * (p+1) + 1:(.5 * (k + 1) * k)])$value)) 

    AL_design <- optim(par = rnorm(batch * p),
        fn = min_eig_graph_batch,
        FI_pre = pre_info_AL, 
        Omega = Omega,beta = beta,
        mu = mu,k = k,p = p,
        control = list(maxit = 3000),method = "L-BFGS-B",lower = -5,upper = 5)
    AL_design <- AL_design$par
    pre_info_AL <- FI_mat_batch(AL_design, 
                                pre_info_AL,
                                Omega,beta,mu,k,p)
    res$AL[i] = min(abs(eigen(
        pre_info_AL[k  * (p+1) + 1:(.5 * (k + 1) * k),
        k * (p+1) + 1:(.5 * (k + 1) * k)])$value))
    cat(i,"\n")

}

par(mfrow = c(1,2))
plot(res$n, res$AL, type = "l",main = paste0("batch=",batch))
lines(res$n, res$random, col = "#1100ff")
lines(res$n, res$rep, col = "#10b30b")

