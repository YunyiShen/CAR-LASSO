# This is a wraper for multi-reg based model, used in horseshoe inference 

# data n*k matrix for data
# design p*k matrix for desing matrix, no intercept
# Bbar: prior mean of regresssion coefficient
# A: prior precision mat for regression coefficnet
# nu d.f. for Sigma
# V: k*k pdf location para for prior on Sigma

CAR_multireg <- function(data,design,n_sample, Bbar=NULL, A = NULL, nu=3, V=NULL){
    k <- ncol(data)
    p <- ncol(design)
    if(is.null(Bbar)) Bbar <- matrix(0,p+1,k)
    if(is.null(A)) A <- diag(1e-8,p+1,p+1)
    if(is.null(V)) V <- 3*diag(2,k,k)

    sample_multireg <- lapply(1:n_sample,
        function(i,data,design,Bbar,A,nu,V){
            temp <- rmultireg(data,cbind(1,design),Bbar,A,nu,V)
            Omega <- solve(temp$Sigma)
            B <- temp$B %*% Omega
            mu <- B[1,]
            B <- (B[-1,])
            return(list(Omega = Omega[upper.tri(Omega,T)],beta = c(B),mu = mu))
        }
        ,data,design,Bbar,A,nu,V)

    Omega_mcmc <- Reduce(rbind,lapply(sample_multireg,function(w){t(w$Omega)}))
    beta_mcmc <- Reduce(rbind,lapply(sample_multireg,function(w){t(w$beta)}))
    mu_mcmc <- Reduce(rbind,lapply(sample_multireg,function(w){t(w$mu)}))
    return(list(Omega = Omega_mcmc,beta = beta_mcmc,mu = mu_mcmc))
}

Multinomial_CAR_multireg <- function(data,design,n_burn_in,n_sample, thin_by, Bbar=NULL, 
                                    A = NULL, nu=3, V=NULL,ns = 1000,m=20,emax=64){
    n <- nrow(data)
    k <- ncol(data)-1
    p <- ncol(design)

    n_store <- floor( n_sample/thin_by )
    i_store <- 1
    if(is.null(Bbar)) Bbar <- matrix(0,p+1,k)
    if(is.null(A)) A <- diag(1e-8,p+1,p+1)
    if(is.null(V)) V <- 3*diag(2,k,k)

    Z_curr <- matrix(rnorm(n*k),n,k)
    mu_curr <- Bbar[1,]
    Omega_curr <- diag(k)
    beta_curr <- t(Bbar[-1,])

    mu_mcmc <- matrix(NA,n_store,k)
    beta_mcmc <- matrix(NA,n_store,k*p)
    Omega_mcmc <- matrix(NA, n_store, floor(.5*k*(k+1)))
    
    # mainloop
    for(i in 1:(n_burn_in + n_sample)){
        #browser()
        temp <- rmultireg(Z_curr,cbind(1,design),Bbar,A,nu,V)
        Omega_curr <- solve(temp$Sigma)
        #browser()
        beta_curr <- temp$B %*% Omega_curr 
        mu_curr <- beta_curr[1,]
        beta_curr <- (beta_curr[-1,])

        update_Z_helper_multinomial_CAR(Z_curr,data,
            design,mu_curr,beta_curr,Omega_curr,k,p,n,ns,m,emax)

        if(i>n_burn_in & ((i-n_burn_in) %% thin_by == 0) ){
            mu_mcmc[i_store,] <- mu_curr
            beta_mcmc[i_store,] <- c(beta_curr)
            Omega_mcmc[i_store,] <- Omega_curr[upper.tri(Omega_curr,T)]

            i_store <- i_store + 1
        }

    }

    return(list(mu = mu_mcmc, beta = beta_mcmc, Omega = Omega_mcmc))

}

Pois_CAR_multireg <- function(data,design,n_burn_in,n_sample, thin_by, Bbar=NULL, 
                                    A = NULL, nu=3, V=NULL,ns = 1000,m=20,emax=64){
    n <- nrow(data)
    k <- ncol(data)
    p <- ncol(design)

    n_store <- floor( n_sample/thin_by )
    i_store <- 1
    if(is.null(Bbar)) Bbar <- matrix(0,p+1,k)
    if(is.null(A)) A <- diag(1e-8,p+1,p+1)
    if(is.null(V)) V <- 3*diag(2,k,k)

    Z_curr <- matrix(rnorm(n*k),n,k)
    mu_curr <- Bbar[1,]
    Omega_curr <- diag(k)
    beta_curr <- t(Bbar[-1,])

    mu_mcmc <- matrix(NA,n_store,k)
    beta_mcmc <- matrix(NA,n_store,k*p)
    Omega_mcmc <- matrix(NA, n_store, floor(.5*k*(k+1)))
    
    # mainloop
    for(i in 1:(n_burn_in + n_sample)){
        #browser()
        temp <- rmultireg(Z_curr,cbind(1,design),Bbar,A,nu,V)
        Omega_curr <- solve(temp$Sigma)
        #browser()
        beta_curr <- temp$B %*% Omega_curr 
        mu_curr <- beta_curr[1,]
        beta_curr <- (beta_curr[-1,])

        update_Z_helper_Pois_CAR(Z_curr,data,
            design,mu_curr,beta_curr,Omega_curr,k,p,n,ns,m,emax)

        if(i>n_burn_in & ((i-n_burn_in) %% thin_by == 0) ){
            mu_mcmc[i_store,] <- mu_curr
            beta_mcmc[i_store,] <- c(beta_curr)
            Omega_mcmc[i_store,] <- Omega_curr[upper.tri(Omega_curr,T)]

            i_store <- i_store + 1
        }

    }

    return(list(mu = mu_mcmc, beta = beta_mcmc, Omega = Omega_mcmc))

}

Pobit_CAR_multireg <- function(data,design,n_burn_in,n_sample, thin_by, Bbar=NULL, 
                                    A = NULL, nu=3, V=NULL){
    n <- nrow(data)
    k <- ncol(data)
    p <- ncol(design)

    n_store <- floor( n_sample/thin_by )
    i_store <- 1
    if(is.null(Bbar)) Bbar <- matrix(0,p+1,k)
    if(is.null(A)) A <- diag(1e-8,p+1,p+1)
    if(is.null(V)) V <- 3*diag(2,k,k)

    Z_curr <- matrix(rnorm(n*k),n,k)
    mu_curr <- Bbar[1,]
    Omega_curr <- diag(k)
    beta_curr <- t(Bbar[-1,])

    mu_mcmc <- matrix(NA,n_store,k)
    beta_mcmc <- matrix(NA,n_store,k*p)
    Omega_mcmc <- matrix(NA, n_store, floor(.5*k*(k+1)))
    
    # mainloop
    for(i in 1:(n_burn_in + n_sample)){
        #browser()
        temp <- rmultireg(Z_curr,cbind(1,design),Bbar,A,nu,V)
        Omega_curr <- solve(temp$Sigma)
        #browser()
        beta_curr <- temp$B %*% Omega_curr 
        mu_curr <- beta_curr[1,]
        beta_curr <- (beta_curr[-1,])

        # update Z, probit model
        update_Z_helper_CAR(Z_curr,data,
            design,mu_curr,beta_curr,Omega_curr,k,p,n,ns,m,emax)

        if(i>n_burn_in & ((i-n_burn_in) %% thin_by == 0) ){
            mu_mcmc[i_store,] <- mu_curr
            beta_mcmc[i_store,] <- c(beta_curr)
            Omega_mcmc[i_store,] <- Omega_curr[upper.tri(Omega_curr,T)]

            i_store <- i_store + 1
        }

    }

    return(list(mu = mu_mcmc, beta = beta_mcmc, Omega = Omega_mcmc))

}