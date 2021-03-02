#' Horseshoe method for graphical structure inference
#' 
#' @details This method fits a linear regression with less informative prior on any parameters and compare the posterior mean with the LASSO result. If LASSO is comparably less than result without sparsity prior, we argue that the edge should be absent
#' @param obj The carlasso_out object from CARlasso
#' @param Bbar Prior mean of regression coefficients, default all 0s
#' @param A Prior precision of regression coefficients, default 1e-8
#' @param nu Prior degree of freedom of the Wishart on precision matrix
#' @param V prior covariance matrix of the Wishart on precision matrix
#' @param thr threshold for horseshoe inference, default 0.5
#' @return A `carlasso_out` object with learned binary adjacency matrix and multi-response linear regression MCMC out put
#' @export
#' @examples 
#' set.seed(42)
#' dt <- simu_AR1()
#' car_res <- CARlasso(y1+y2+y3+y4+y5~x1+x2+x3+x4+x5, data = dt, adaptive = TRUE)
#' car_res <- horseshoe(car_res)
#' plot(car_res)


horseshoe <- function(obj, Bbar=NULL, A = NULL, nu=3, V=NULL, thr = 0.5 ){
    y <- obj$data$response
    design <- obj$data$design
    ns <- obj$settings$ns
    m <- obj$settings$m
    emax <- obj$settings$emax
    if(obj$settings$link == "identity"){
        multireg_res <- CAR_multireg(y,design,nrow(obj$MCMC_output$beta),
                                     Bbar, A, nu, V)
    }

    if(obj$settings$link == "probit"){
        multireg_res <- Pobit_CAR_multireg(y,design,obj$settings$n_burn_in,obj$settings$n_iter,
                                     obj$settings$thin_by, 
                                     Bbar, A, nu, V)
    }

    if(obj$settings$link == "log"){
        multireg_res <- Pois_CAR_multireg(y,design,obj$settings$n_burn_in,obj$settings$n_iter,
                                     obj$settings$thin_by, 
                                     Bbar, A, nu, V, ns, m, emax)
    }

    if(obj$settings$link == "logit"){
        multireg_res <- Multinomial_CAR_multireg(y,design,obj$settings$n_burn_in,obj$settings$n_iter,
                                     obj$settings$thin_by, 
                                     Bbar, A, nu, V, ns, m, emax)
    }


    graph_multireg <- get_graph(multireg_res, k = nrow(obj$point_est$Omega))
    B_multireg <- matrix(colMeans(multireg_res$beta),nrow = nrow(obj$point_est$beta))

    horseshoe_binary <- list(Omega_binary = abs(obj$point_est$Omega/graph_multireg)>thr,
                            B_binary = abs(abs(obj$point_est$beta/B_multireg)>thr))

    obj$horseshoe_binary <- horseshoe_binary
    obj$multireg_mcmc <- multireg_res
    return(obj)

}


