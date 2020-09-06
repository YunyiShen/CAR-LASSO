g_model1 <- function(k, rho=.7){
    temp <- matrix(rep(1:k,k),ncol = k)
    Sigma <- rho ^ (abs(temp-t(temp)))
    Omega <- solve(Sigma)
    Omega <- Omega * (abs(Omega)>1e-15)
    return(list(Sigma = Sigma, Omega = Omega))
}


g_model2 <- function(k, rho = .5){
    temp <- matrix(rep(1:k,k),ncol = k)
    Omega <- rho ^ (abs(temp-t(temp))) * (abs(temp-t(temp)) <= 2)
    return(list(Omega = Omega, Sigma = solve(Omega)))
}

g_model3 <- function(k, rho = .5){
    row_ind <- matrix(rep(1:k,k),ncol = k)
    col_ind <- t(row_ind)

    Sigma <- diag(1,k,k)

    Sigma[(1<=row_ind ) &
          (1<=col_ind ) &
          (row_ind != col_ind) & 
          (row_ind <= k/2) & 
          (col_ind <= k/2)] <- rho
    Sigma[((1+k/2)<=row_ind ) & 
          ((1+k/2)<=col_ind ) &
          (row_ind != col_ind) & 
          (row_ind <= 10) & 
          (col_ind <= 10)] <- rho
    Omega <- solve(Sigma)
    
    return(list(Sigma = Sigma, Omega = Omega))
}

g_model4 <- function(k, rho=.1){
    
    Omega <- diag(1,k,k)
    Omega[2:k,1] <- rho
    Omega[1,2:k] <- rho
    Sigma <- solve(Omega)
    return(list(Sigma = Sigma, Omega = Omega))
}

g_model5 <- function(k, rhos = c(2,1,.9)){
    row_ind <- matrix(rep(1:k,k),ncol = k)
    col_ind <- t(row_ind)
    Omega <- diag(rhos[1],k,k)
    Omega[abs(row_ind-col_ind)==1] <- rhos[2]
    Omega[1,k] <- Omega[k,1] <- rhos[3]
    return(list(Sigma = solve(Omega), Omega = Omega))
}

g_model6 <- function(k, rhos = c(2,1)){
    Omega <- matrix(rhos[2], k, k)
    diag(Omega) <- rhos[1]

    return(list(Sigma = solve(Omega),Omega = Omega))
}