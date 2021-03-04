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

    res <- CAR_multireg_cpp(data,design,n_sample, Bbar, A, nu, V)
    return(res)
}

Multinomial_CAR_multireg <- function(data,design,n_burn_in,n_sample, thin_by, Bbar=NULL, 
                                    A = NULL, nu=3, V=NULL,ns = 1000,m=20,emax=64){
    n <- nrow(data)
    k <- ncol(data)-1
    p <- ncol(design)

    
    if(is.null(Bbar)) Bbar <- matrix(0,p+1,k)
    if(is.null(A)) A <- diag(1e-8,p+1,p+1)
    if(is.null(V)) V <- 3*diag(2,k,k)

    res <- Multinomial_CAR_multireg_cpp(data,design,n_burn_in,n_sample, thin_by, Bbar, 
                                    A, nu, V,ns,m,emax)
    return(res)

}

Pois_CAR_multireg <- function(data,design,n_burn_in,n_sample, thin_by, Bbar=NULL, 
                                    A = NULL, nu=3, V=NULL,ns = 1000,m=20,emax=64){
    n <- nrow(data)
    k <- ncol(data)
    p <- ncol(design)

    
    if(is.null(Bbar)) Bbar <- matrix(0,p+1,k)
    if(is.null(A)) A <- diag(1e-8,p+1,p+1)
    if(is.null(V)) V <- 3*diag(2,k,k)

    res <- Pois_CAR_multireg_cpp(data,design,n_burn_in,n_sample, thin_by, Bbar, 
                                    A, nu, V,ns,m,emax)
    return(res)

}

Probit_CAR_multireg <- function(data,design,n_burn_in,n_sample, thin_by, Bbar=NULL, 
                                    A = NULL, nu=3, V=NULL){
    n <- nrow(data)
    k <- ncol(data)
    p <- ncol(design)

    if(is.null(Bbar)) Bbar <- matrix(0,p+1,k)
    if(is.null(A)) A <- diag(1e-8,p+1,p+1)
    if(is.null(V)) V <- 3*diag(2,k,k)

    res <- Probit_CAR_multireg_cpp(data,design,n_burn_in,n_sample, thin_by, Bbar, 
                                    A, nu, V)

    return(res)

}
