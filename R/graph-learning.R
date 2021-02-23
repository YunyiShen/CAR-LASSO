


horseshoe.carlasso_out(obj, Bbar=NULL, A = NULL, nu=3, V=NULL ){
    y <- obj$data$response
    design <- obj$data$design
    if(obj$settings$link = "identity"){
        multireg_res <- CAR_multireg(y,design,nrow(obj$MCMC_output$beta),
                                     Bbar, A, nu, V)
    }

    if(obj$settings$link = "probit"){
        multireg_res <- Pobit_CAR_multireg(y,design,obj$settings$n_burn_in,obj$settings$n_iter,
                                     obj$settings$thin_by, 
                                     Bbar, A, nu, V)
    }




}


