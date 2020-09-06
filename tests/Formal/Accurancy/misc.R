get_graph <- function(CAR_sample,k){
  Omega <- matrix(0,k,k)
  Omega[upper.tri(Omega,T)] = apply(CAR_sample$Omega,2,mean)
  Omega <- Omega+t(Omega)
  diag(Omega) <- 0.5 * diag(Omega)
  return(Omega)
}


get_beta <- function(CAR_sample,p){
    matrix(colMeans(CAR_sample$beta),nrow = p)
}

multireg_Sample <- function(data, Design, k, p, n_sample = 10000){
    sample_multireg <- lapply(1:n_sample,function(i) 
        rmultireg(data,cbind(1,Design),
            matrix(0,p+1,k),
            diag(1e-8,p+1,p+1),
            3,3*diag(2,k,k))
        )
    Beta_multireg <- lapply(sample_multireg,function(w) t(solve(w$Sigma,t(w$B[-1,]))))
    Omega_multireg <- lapply(sample_multireg,function(w) w$Sigma)
    Omega_multireg <- lapply(Omega_multireg,solve)

    Beta_hat <- Reduce("+",Beta_multireg)/n_sample
    multireg_Graph <- Reduce("+",Omega_multireg)/n_sample
    return(list(beta = Beta_hat, Omega = multireg_Graph))
}

makedata <- function(Sigma,design,beta,n){
    Xbeta <- design %*% beta
    Z <- matrix(0,nrow = n, ncol = nrow(Sigma))
    for(i in 1:n){
        Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]),Sigma)
    }
    return(Z)
}

get_counts_Omega <- function(graph_hat,graph){
    graph <- graph!=0
    TP <- sum(graph_hat[upper.tri(graph_hat)] & graph[upper.tri(graph)])
    TN <- sum(!graph_hat[upper.tri(!graph_hat)] & !graph[upper.tri(!graph)])
    FP <- sum(graph_hat[upper.tri(graph_hat)] & !graph[upper.tri(!graph)])
    FN <- sum(!graph_hat[upper.tri(!graph_hat)] & graph[upper.tri(graph)])
    return(c(TP,TN,FP,FN))
}

get_counts_beta <- function(beta_hat,beta){
    beta <- beta != 0
    TP <- sum(beta_hat & beta)
    TN <- sum(!beta_hat & !beta)
    FP <- sum(beta_hat & !beta)
    FN <- sum(!beta_hat & beta)
    return(c(TP,TN,FP,FN))
}


log_l2 <- function(a,b){
    log(sum((a-b)^2))
}