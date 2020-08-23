library(Rcpp)
library(RcppArmadillo)
library(Matrix)
sourceCpp("./src/CAR_FI_helper.cpp")

min_eig_graph <- function(design,FI_pre,Omega,beta,mu,k,p){
    design <- t(design)
    temp <- FI_pre + CAR_FI(design,Omega,beta ,mu,k,p)
    return(-log(min(abs(eigen(
            temp[(k ) * (p+1) + 1:(.5 * (k + 1) * k),
                 (k ) * (p+1) + 1:(.5 * (k + 1) * k)]
        )$value))))
}

det_graph <- function(design,FI_pre,Omega,beta,mu,k,p){
    design <- t(design)
    temp <- FI_pre + CAR_FI(design,Omega,beta ,mu,k,p)
    return(-log(det(
            temp[(k ) * (p+1) + 1:(.5 * (k + 1) * k),
                 (k ) * (p+1) + 1:(.5 * (k + 1) * k)]
        )))
}

min_eig_all <- function(design,FI_pre,Omega,beta,mu,k,p){
    design <- t(design)
    temp <- FI_pre + CAR_FI(design,Omega,beta ,mu,k,p)
    return(-min(abs(eigen(
            temp
        )$value)))
}

k <- 10
p <- 10

beta <- matrix(rnorm(k * p) , p,k)
mu <- rnorm(k)

B <- rsparsematrix(k,k,0.2)
omega <- diag(rgamma(k,3,.1))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
#Omega <- omega
#diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)


init_design <- matrix(rnorm(p*k),k,p)
init_info <- CAR_FI(t(init_design[1,]),Omega,beta ,mu,k,p)
for(i in 2:k){
    init_info <- CAR_FI(t(init_design[i,]),Omega,beta ,mu,k,p) + init_info
}
pre_info_AL <- init_info
pre_info_rep <- init_info
pre_info_rand <- init_info

AL_round <- 100
res <- data.frame(n = 1:AL_round,rep = NA, random = NA, AL = NA)

for(i in 1:AL_round){
    pre_info_rep <- pre_info_rep + 
        CAR_FI(t(init_design[1,]),Omega,beta ,mu,k,p ) 
    res$rep[i] = log(det(pre_info_rep[(k ) * (p+1) + 1:(.5 * (k + 1) * k),(k ) * (p+1) + 1:(.5 * (k + 1) * k)]))

    pre_info_rand <- pre_info_rand + 
        CAR_FI(t(rnorm(p)),Omega,beta ,mu,k,p)
    res$random[i] = log(det(pre_info_rand[(k ) * (p+1) + 1:(.5 * (k + 1) * k),(k ) * (p+1) + 1:(.5 * (k + 1) * k)]))

    AL_design <- optim(par = init_design[1,],
        fn = det_graph,
        FI_pre = pre_info_AL, 
        Omega = Omega,beta = beta,mu = mu,k = k,p = p)
    AL_design <- t(AL_design$par)
    pre_info_AL <- pre_info_AL + 
        CAR_FI(AL_design,Omega,beta ,mu,k,p)
    res$AL[i] = log(det(pre_info_AL[(k ) * (p+1) + 1:(.5 * (k + 1) * k),(k ) * (p+1) + 1:(.5 * (k + 1) * k)]))
    cat(i,"\n")

}

plot(res$n, res$AL, type = "l")
lines(res$n, res$random, col = "red")
lines(res$n, res$rep, col = "blue")