# This file calculate information for optimal design 


# Bayesian D-optimal
al_info_D <- function(design_vec,CAR_sample,k,pp,
                           nrep, pan_scale=.1,pan_amp=.01, 
                           summary = "mean"){
    p <- pp
    # point estimation 
    Omega_hat <- get_graph(CAR_sample,k, summary)
    beta_hat <- matrix(apply(CAR_sample$beta,2,summary),nrow = p)
    mu_hat <- apply(CAR_sample$mu,2,summary)

    # Laplace approximation of the prior, i.e. last step's posterior 
    Emp_cov_prior <- cov(Reduce(cbind,CAR_sample[1:3]))

    # information: basically percision of prior/posterior 
    Emp_info_prior <- solve(Emp_cov_prior)

    # get design matrix, due to repeat
    designmat <- matrix(design_vec,ncol = p)
    designmat <- lapply(1:nrep,function(i) designmat)
    designmat <- Reduce(rbind,designmat)

    # calculate information matrix 
    infomat <- CAR_FI_para(designmat,Omega_hat,beta_hat,mu_hat,k,p=p) + 
            Emp_info_prior
    return(-log(det(infomat))+ 
        sum(pan_amp * exp(pan_scale * designmat^2))/length(designmat)) # add penalty 
}

al_info_A <- function(design_vec,CAR_sample,k,
                          pp,nrep,pan_scale=.1,pan_amp=.01,
                          summary = "mean"){
        p <- pp
    # point estimation 
    Omega_hat <- get_graph(CAR_sample,k, summary)
    beta_hat <- matrix(apply(CAR_sample$beta,2,summary),nrow = p)
    mu_hat <- apply(CAR_sample$mu,2,summary)

    # Laplace approximation of the prior, i.e. last step's posterior 
    Emp_cov_prior <- cov(Reduce(cbind,CAR_sample[1:3]))

    # information: basically percision of prior/posterior 
    Emp_info_prior <- solve(Emp_cov_prior)



    designmat <- matrix(design_vec,ncol = p)
    designmat <- lapply(1:nrep,function(i) designmat)
    designmat <- Reduce(rbind,designmat)
  
    infomat <- CAR_FI_para(designmat,Omega_hat,beta_hat,mu_hat,k,p) + 
        Emp_info_prior
    return(log(sum(diag(solve(infomat)))) + 
           sum(pan_amp * exp(pan_scale * designmat^2)) / length(designmat))
}