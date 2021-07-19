#' Gibbs sampler for Conditional Autoregressive LASSO and extensions
#' 
#' @description Main sampling algorithm of CAR-LASSO model
#'
#' @param formula A double sided formula with response at left hand side and predictors at right hand side
#' @param data A data.frame with all response and predictors, row as observations
#' @param link String name of link function? Currently can be "identity" for normal response, "probit" for binary, "log" for counting, "logit" for compositional. Note that when use "logit", the last response will be used as reference. 
#' @param adaptive Bool, whether run the adaptive version of the model
#' @param r_beta Hyper-parameter for regression coefficient, shape parameter of Gamma, if adaptive, should have row number same as number of predictors while column number of responses
#' @param delta_beta Hyper-parameter for regression coefficient, rate parameter of Gamma, if adaptive, should have row number same as number of predictors while column number of responses
#' @param r_Omega Hyper-parameter for precision matrix, shape parameter of Gamma. If adaptive, can be a matrix with same size as precision matrix, if this is the case, only upper triangular part without diagonal will be used, or can be a vector whose size was the upper triangular part of precision matrix, if non-adaptive, a scalar.
#' @param delta_Omega Hyper-parameter for precision matrix, rate parameter of Gamma, If adaptive, can be a matrix with same size as precision matrix, if this is the case, only upper triangular part without diagonal will be used, or can be a vector whose size was the upper triangular part of precision matrix, if non-adaptive, a scalar.
#' @param lambda_diag adaptive only hyper-parameter for penalties on diagonal entries of Omega, should have dimension k and non-negative
#' @param n_iter Number of sampling iterations (i.e. after burn in) for the Gibbs sampler
#' @param n_burn_in Number of burn in iterations for the Gibbs sampler
#' @param thin_by Final sample was thin by this number
#' @param ns parameter for ARS, maximum number of hulls, only used when link is "log" and "logit"
#' @param m parameter for ARS, initial number of hulls, only used when link is "log" and "logit"
#' @param emax parameter for ARS, tolerance for small values being 0, larger meaning we tolerate smaller values, only used when link is "log" and "logit"
#' @param progress Bool, whether report progress from C++
#' @param verbos Bool, whether show warnings and messages.
#' 
#' @return A `carlasso_out` object with elements: 
#' \itemize{
#'    \item{`$point_est`}{
#'        \itemize{
#'          \item{`$Omega`}{: Posterior mean of precision matrix}
#'          \item{`$beta`}{: Posterior mean of regression coefficient}
#'          \item{`$CAR`}{
#'            \itemize{
#'              \item{`$C`}{: The conditional regression coefficients among responses}
#'              \item{`$B`}{: The conditional regression coefficients between response and predictors}
#'              \item{`$M`}{: The conditional variance}
#'            }
#'          }
#'        }
#'    }
#'    \item{`$nodes`}{
#'        \itemize{
#'            \item{`$responses`}{: node name of responses}
#'            \item{`$predictors`}{: node name of predictors}
#'        }
#'    }
#' 
#'    \item{`$data`}{
#'        \itemize{
#'            \item{`$response`}{: response matrix}
#'            \item{`$design`}{: design matrix}
#'        }
#'    }
#' 
#'    \item{`$settings`}{: all settings sent to the algorithm, exclude data}
#'    \item{`$MCMC_output`}{
#'        \itemize{
#'            \item{`$beta`}{: A coda::mcmc object, each row was an MCMC sample of the (column) vectorization of regression coefficient B}
#'            \item{`$mu`}{: A coda::mcmc object, each row was an MCMC sample of the mean vector}
#'            \item{`$Omega`}{: A coda::mcmc object, each row was an MCMC sample of the upper triangular part (with diagonal) of precision matrix Omega}
#'            \item{`$lambda`}{: \strong{Non-adaptive only}, A coda::mcmc object, first column was the shrinkage parameter lambda for regression coefficient and the second column was shrinkage parameter lambda for precision matrix}
#'            \item{`$lambda_beta`}{: \strong{Adaptive only}, A coda::mcmc object, each row was an MCMC sample of the (column) vectorization of shrinkage parameter for regression coefficient B}
#'            \item{`$lambda_Omega`}{: \strong{Adaptive only}, A coda::mcmc object, each row was an MCMC sample of the shrinkage parameter for the upper triangular part (without diagonal) of precision matrix Omega}
#'        }
#'    }
#' }
#' 
#' 
#' 
#' @examples
#' set.seed(42)
#' dt <- simu_AR1()
#' car_res <- CARlasso(y1+y2+y3+y4+y5~x1+x2+x3+x4+x5, data = dt, adaptive = TRUE)
#' plot(car_res,tol = 0.05)
#' # with horseshoe inference
#' car_res <- horseshoe(car_res)
#' plot(car_res)
#' 
#' 


CARlasso <- function(formula, # a double sided formula needed, e.g. x+y~a+b
                     data, link = "identity",
                     adaptive = FALSE,
                     r_beta = ifelse(adaptive,0.01,1), 
                     delta_beta = ifelse(adaptive,1e-6,0.01),
                     r_Omega = ifelse(adaptive,0.01,1), 
                     delta_Omega = ifelse(adaptive,1e-6,0.01),
                     lambda_diag = 0,
                     n_iter = 2000,
                     n_burn_in = 1000, thin_by = 10,
                     ns = 1000, m=20, emax=64, 
                     progress = TRUE, verbos = TRUE) {
  # some warning messages
  err_no_predictor <- "No predictor supplied, do you want to try bGlasso?.\n\n"
  warr_centering <- "Predictors will be centered.\n\n"

  # check links
  if (!(link %in% c("identity", "probit", "log", "logit"))) {
    stop("Currently only implemented identity (normal), 
       log (Poisson) logit (multi-nomial) and probit (bernoulli)")
  }

  # omit NAs
  if (!all(!is.na(data)) & verbos) {
    warning("NAs in data are omitted")
    data <- na.omit(data)
  }

  # get matrices from data
  design <- model.matrix(formula, data)
  # design matrix, no intercept:
  design <- design[, colnames(design) != "(Intercept)"]
  design <- as.matrix(design)


  response <- formula
  response[[3]] <- formula[[2]]
  y <- model.matrix(response, data) # response matrix.
  y <- y[, colnames(y) != "(Intercept)"]
  rm(response)

  # center predictors
  if (ncol(design) == 0) {
    stop(err_no_predictor)
  }

  if (verbos) {
    cat(warr_centering)
  }

  design <- apply(design, 2, function(w) {
    (w - mean(w)) / sd(w)
  })

  # get dimensions
  n <- nrow(design)
  p <- ncol(design)
  k <- ncol(y) - (link == "logit") # careful for multinomial response

  ### check dimension of r_beta
  dimname1 <- c(
    "nrow r_beta", "ncol r_beta",
    "nrow delta_beta", "ncol delta_beta"
  )
  if (!all(c(c(r_beta), c(delta_beta)) > 0)) {
    stop("Hyperparameters for beta shrinkage must be positive\n\n")
  }
  if (!adaptive) {
    if (is.null(r_beta)) r_beta <- 1
    if (is.null(delta_beta)) delta_beta <- 0.01
    if (verbos & (length(r_beta) > 1 | length(delta_beta) > 1)) {
      cat("Algorithm set to be non-adapive, will take the first entry of hyperprior for beta shrinkage\n\n")
    }
    r_beta <- r_beta[1]
    delta_beta <- delta_beta[1]
  }
  else {
    if (is.null(r_beta)) r_beta <- 0.01
    if (is.null(delta_beta)) delta_beta <- 1e-6
    if ((length(r_beta) == 1 & length(delta_beta) == 1)) {
      if (verbos) cat("Algorithm set to be adapive. Assuming all hyper parameters are the same for beta \n\n")
      r_beta <- matrix(r_beta, p, k)
      delta_beta <- matrix(delta_beta, p, k)
    }
    else {
      dims <- c(nrow(r_beta), ncol(r_beta), nrow(delta_beta), ncol(delta_beta))
      prop_dim <- c(p, k, p, k)
      mismatch <- dimname1[dims != prop_dim]
      if (length(mismatch) > 0) {
        errmsg <- paste(
          "Dimension mismatch for hyper prior of beta shrinkage: ",
          paste(mismatch, collapse = " "), "\n\n"
        )
        stop(errmsg)
      }
    }
  }
  ## end checking r_beta

  ## checking r_Omega, delta_Omega
  if (is.matrix(r_Omega)) {
    if (verbos) cat("Supplied matrix for hyper parameter for r_Omega, will only take upper triangular part\n\n")
    r_Omega <- c(r_Omega[upper.tri(r_Omega)])
  }

  if (is.matrix(delta_Omega)) {
    if (verbos) cat("Supplied matrix for hyper parameter for delta_Omega, will only take upper triangular part\n\n")
    delta_Omega <- c(delta_Omega[upper.tri(delta_Omega)])
  }

  if (!adaptive) {
    if (is.null(r_Omega)) r_Omega <- 1
    if (is.null(delta_Omega)) delta_Omega <- 0.01
    if (verbos & (length(r_Omega) > 1 | length(delta_Omega) > 1)) {
      cat("Algorithm set to be non-adapive, will take the first entry of hyper prior for Omega shrinkage\n\n")
    }
    r_Omega <- r_Omega[1]
    delta_Omega <- delta_Omega[1]
  }
  else {
    if (is.null(r_Omega)) r_Omega <- 0.01
    if (is.null(delta_Omega)) delta_Omega <- 1e-6
    if (is.null(lambda_diag)) lambda_diag <- 0
    if ((length(r_Omega) == 1 & length(delta_Omega) == 1)) {
      if (verbos) cat("Algorithm set to be adapive. Assuming all hyper parameters are the same for Omega's off diagonal entries \n\n")
      r_Omega <- rep(r_Omega, .5 * (k - 1) * k)
      delta_Omega <- rep(delta_Omega, .5 * (k - 1) * k)
    }
    else {
      if (length(r_Omega) != .5 * (k - 1) * k |
        length(delta_Omega) != .5 * (k - 1) * k) {
        errmsg <- "Dimension mismatch for hyper prior of Omega shrinkage \n\n"
        stop(errmsg)
      }
    }
    if ((length(lambda_diag) == 1 )) {
      if (verbos) cat("Algorithm set to be adaptive. Assuming priors are all the same for Omega's diagonals \n\n")
      lambda_diag <- rep(lambda_diag, k)
    }
    else {
       if(length(lambda_diag) != k ){
         errmsg <- "Dimension mismatch for hyper prior of Omega diagonal shrinkage \n\n"
        stop(errmsg)
       }
    }
  }
  ## end checking Omega hyper parameters

  
  ## Main algorithms

  if (link == "identity") {
    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    if (adaptive) {
      res <- CAR_ALASSO_Cpp(
        y, design,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, lambda_diag,
        progress
      )
    }
    else {
      res <- CAR_LASSO_Cpp(
        y, design,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, progress
      )
    }
  }

  if (link == "log") {
    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    if (adaptive) {
      res <- CAR_ALASSO_hir_Cpp(
        y, design, link = 1,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, lambda_diag,
        ns, m, emax, 
        progress
      )
    }
    else {
      res <- CAR_LASSO_hir_Cpp(
        y, design, link = 1,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega,
        ns, m, emax, 
        progress
      )
    }
  }

  if (link == "logit") {
    if (verbos) cat("Last response will be used as reference group\n\n")
    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    if (adaptive) {
      res <- CAR_ALASSO_hir_Cpp(
        y, design, link = 2,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, lambda_diag,
        ns, m, emax, 
        progress
      )
    }
    else {
      res <- CAR_LASSO_hir_Cpp(
        y, design, link = 2,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega,
        ns, m, emax, 
        progress
      )
    }
  }
  if (link == "probit") {
    unique_values <- apply(y, 2, function(w) {
      length(unique(w))
    })
    if (!all(unique_values == 2)) {
      stop("Response has multiple unique values, cannot use probit link, do you mean logit?\n")
    }

    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    if (adaptive) {
      res <- CAR_ALASSO_hir_Cpp(
        y, design, link = 3,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, lambda_diag,
        ns, m, emax, 
        progress
      )
    }
    else {
      res <- CAR_LASSO_hir_Cpp(
        y, design, link = 3,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega,
        ns, m, emax, 
        progress
      )
    }
  }

  omega_post <- get_graph(res,k)
  b_post <- matrix(colMeans(res$beta),p,k)
  CAR_post <- get_CAR_MB(b_post,omega_post)

  point_est <- list(Omega = omega_post, beta = b_post, CAR = CAR_post)

  res <- lapply(res, coda::mcmc)

  settings <- list(formula = formula, link = link, adaptive = adaptive,
                     r_beta = r_beta , delta_beta = delta_beta , r_Omega = r_Omega, 
                     delta_Omega = delta_Omega, lambda_diag = lambda_diag, n_iter = n_iter,
                     n_burn_in = n_burn_in, thin_by = thin_by, ns=ns,m=m,emax=emax,progress = progress, verbos = verbos)

  nodes <- list(response = colnames(y), predictors = colnames(design))
  res <- list(point_est = point_est, nodes = nodes,  
              data = list(response = y, design = design),
              settings = settings, MCMC_output = res)

  class(res) <- "carlasso_out"

  if(verbos) cat("\ndone\n\n")

  return(res)
}

#' Gibbs sampler for Bayesian Graphical LASSO and extensions
#' 
#' @description Main sampling algorithm of Glasso model, note that the mean is in CAR parameterization
#'
#' @param data A data.frame with all response, row as observations
#' @param link String name of link function? Currently can be "identity" for normal response, "probit" for binary, "log" for counting, "logit" for compositional. Note that when use "logit", the last response will be used as reference. 
#' @param r_Omega Hyper-parameter for precision matrix, shape parameter of Gamma. Should be a scalar 
#' @param delta_Omega Hyper-parameter for precision matrix, rate parameter of Gamma. Shoule be a scalar
#' @param n_iter Number of sampling iterations (i.e. after burn in) for the Gibbs sampler
#' @param n_burn_in Number of burn in iterations for the Gibbs sampler
#' @param thin_by Final sample was thin by this number
#' @param ns parameter for ARS, maximum number of hulls, only used when link is "log" and "logit"
#' @param m parameter for ARS, initial number of hulls, only used when link is "log" and "logit"
#' @param emax parameter for ARS, tolerance for small values being 0, larger meaning we tolerate smaller values, only used when link is "log" and "logit"
#' @param progress Bool, whether report progress from C++
#' @param verbos Bool, whether show warnings and messages.
#' 
#' @return A `bglasso_out` object with elements: 
#' \itemize{
#'    \item{`$point_est`}{
#'        \itemize{
#'          \item{`$Omega`}{: Posterior mean of precision matrix}
#'        }
#'    }
#'    \item{`$nodes`}{
#'        \itemize{
#'            \item{`$responses`}{: node name of responses}
#'        }
#'    }
#' 
#'    \item{`$data`}{
#'        \itemize{
#'            \item{`$response`}{: response matrix}
#'        }
#'    }
#' 
#'    \item{`$settings`}{: all settings sent to the algorithm, exclude data}
#'    \item{`$MCMC_output`}{
#'        \itemize{
#'            \item{`$mu`}{: A coda::mcmc object, each row was an MCMC sample of the mean vector}
#'            \item{`$Omega`}{: A coda::mcmc object, each row was an MCMC sample of the upper triangular part (with diagonal) of precision matrix Omega}
#'            \item{`$lambda`}{:  A coda::mcmc object, first column was the shrinkage parameter lambda for regression coefficient and the second column was shrinkage parameter lambda for precision matrix}
#'        }
#'    }
#' }
#' 
#' 
#' 
#' @examples
#' set.seed(42)
#' dt <- simu_AR1()
#' glassores <- bGlasso(data = dt[,1:5])
#' plot(glassores) 
#' 


bGlasso <- function(data, link = "identity",
                     r_Omega = 1, 
                     delta_Omega = 0.01,
                     n_iter = 2000,
                     n_burn_in = 1000, thin_by = 10,
                     ns = 1000, m=20, emax=64, 
                     progress = TRUE, verbos = TRUE) {

  # check links
  if (!(link %in% c("identity", "probit", "log", "logit"))) {
    stop("Currently only implemented identity (normal), 
       log (Poisson) logit (multi-nomial) and probit (bernoulli)")
  }

  # omit NAs
  if (!all(!is.na(data)) & verbos) {
    warning("NAs in data are omitted")
    data <- na.omit(data)
  }




  

  y <- as.matrix(data)
  # get dimensions
  n <- nrow(y)
  k <- ncol(y) - (link == "logit") # careful for multinomial response


  
    if (is.null(r_Omega)) r_Omega <- 1
    if (is.null(delta_Omega)) delta_Omega <- 0.01
    if (verbos & (length(r_Omega) > 1 | length(delta_Omega) > 1)) {
      cat("Multiple hyper parameters supplied, will take the first entry of hyper prior for Omega shrinkage\n\n")
    }
    r_Omega <- r_Omega[1]
    delta_Omega <- delta_Omega[1]
  
  ## end checking Omega hyper parameters

  
  ## Main algorithms

  if (link == "identity") {
    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    
      res <- Intercept_Graphical_LASSO_Cpp(
        y,n_iter, n_burn_in,
        thin_by,
        r_Omega, delta_Omega,
        progress
      )
  }

  if (link == "log") {
    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    
    res <- Intercept_Graphical_LASSO_hir_Cpp(
        y, 1,
        n_iter, n_burn_in,
        thin_by,
        r_Omega, delta_Omega,
        ns, m, emax, 
        progress
    )
  }

  if (link == "logit") {
    if (verbos) cat("Last response will be used as reference group\n\n")
    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    res <- Intercept_Graphical_LASSO_hir_Cpp(
        y, 2,
        n_iter, n_burn_in,
        thin_by,
        r_Omega, delta_Omega,
        ns, m, emax, 
        progress
    )
    
  }
  if (link == "probit") {
    unique_values <- apply(y, 2, function(w) {
      length(unique(w))
    })
    if (!all(unique_values == 2)) {
      stop("Response has multiple unique values, cannot use probit link, do you mean logit?\n")
    }

    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    
    
    res <- Intercept_Graphical_LASSO_hir_Cpp(
        y, 3,
        n_iter, n_burn_in,
        thin_by,
        r_Omega, delta_Omega,
        ns, m, emax, 
        progress
    )
  }

  omega_post <- get_graph(res,k)

  point_est <- list(Omega = omega_post)

  res <- lapply(res, coda::mcmc)

  settings <- list(link = link,
                    r_Omega = r_Omega, 
                     delta_Omega = delta_Omega, n_iter = n_iter,
                     n_burn_in = n_burn_in, thin_by = thin_by, ns=ns,m=m,emax=emax,progress = progress, verbos = verbos)

  nodes <- list(response = colnames(y))
  res <- list(point_est = point_est, nodes = nodes,  
              data = list(response = y),
              settings = settings, MCMC_output = res)

  class(res) <- "bglasso_out"

  if(verbos) cat("\ndone\n\n")

  return(res)
}


