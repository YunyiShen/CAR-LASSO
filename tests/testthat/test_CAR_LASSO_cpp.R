test_that("CARlasso.cpp",{
    library(Matrix)
    library(RcppArmadillo)
    library(Rcpp)
    library(RcppProgress)


    k <- 11
    n <- 8000
    p <- 1
    B <- as.matrix( Matrix::rsparsematrix(k,k,0.1))
    Omega <- diag(rgamma(k,5,.1))
    IminusB<- 0-B
    diag(IminusB) <- 1+diag(IminusB)
    Omega <- t(IminusB) %*% Omega %*% IminusB

    #Omega <- as.matrix(Omega)
    

    Sigma <- solve(Omega)

    Design <- matrix(rnorm(n*p,0,1),n,p)


    beta <- matrix(rnorm(p*k,1,1),p,k)

    mu <-  rnorm(k)

    Xbeta <- Design %*% beta

    Z <- matrix(NA,n,k)

    for( i in 1:n ){
        Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+mu),Sigma)
    }
    CAR_LASSO_Cpp(Z,  Design, n_iter = 50, 
                          n_burn_in = 50, thin_by = 1, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T)

})