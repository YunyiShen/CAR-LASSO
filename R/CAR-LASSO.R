#' Gibbs sampler for Conditional Auto Regression LASSO and extensions
#'
#' @param formula A double sided formula with response at left hand side and predictors at right hand side
#' @param data A data.frame with all response and predictors, row as observations
#' @param link String name of link function? Currently can be "identity" for normal response, "probit" for binary, "log" for counting, "logit" for compositional. Note that when use "logit", the last response will be used as reference. 
#' @param adaptive Bool, whether run the adaptive version of the model
#' @param r_beta Hyper-parameter for regression coefficient, shape parameter of Gamma, if adaptive, should have row number same as number pf predictors while column number of responses
#' @param delta_beta Hyper-parameter for regression coefficient, rate parameter of Gamma, if adaptive, should have row number same as number pf predictors while column number of responses
#' @param r_Omega Hyper-parameter for precision matrix, shape parameter of Gamma. If adaptive, can be a matrix with same size as precision matrix, if this is the case, only upper triangular part without diagonal will be used, or can be a vector whose size was the upper triangular part of precision matrix, if non-adaptive, a scalar.
#' @param delta_Omega Hyper-parameter for precision matrix, rate parameter of Gamma, If adaptive, can be a matrix with same size as precision matrix, if this is the case, only upper triangular part without diagonal will be used, or can be a vector whose size was the upper triangular part of precision matrix, if non-adaptive, a scalar.
#' 
#' @param n_iter Number of sampling iterations (i.e. after burn in) for the Gibbs sampler
#' @param n_burn_in Number of burn in iterations for the Gibbs sampler
#' @param thin_by Final sample was thin by this number
#' @param prograss Bool, whether report progress from C++
#' @param verbos Bool, whether show warnings and messages.
#' 
#' @return A list of MCMC samples. 
#' \item{beta}{A coda::mcmc object, each row was an MCMC sample of the (column) vectorization of regression coefficient B}
#' \item{mu}{A coda::mcmc object, each row was an MCMC sample of the mean vector}
#' \item{Omega}{A coda::mcmc object, each row was an MCMC sample of the upper triangular part (with diagonal) of precision matrix Omega}
#' \item{lambda}{\strong{Non-adaptive only}, A coda::mcmc object, first column was the shrinkage parameter lambda for regression coefficient and the second column was shrinkage parameter lambda for precision matrix}
#' \item{lambda_beta}{\strong{Adaptive only}, A coda::mcmc object, each row was an MCMC sample of the (column) vectorization of shrinkage parameter for regression coefficient B}
#' \item{lambda_Omega}{\strong{Adaptive only}, A coda::mcmc object, each row was an MCMC sample of the shrinage parameter for the upper triangular part (without diagonal) of precision matrix Omega}
#' 
#' @examples
#' 
#' 



CARlasso <- function(formula, # a double sided formula needed, e.g. x+y~a+b
                     data, link = "identity",
                     adaptive = F,
                     r_beta = ifelse(adaptive,0.01,1), 
                     delta_beta = ifelse(adaptive,1e-6,0.01),
                     r_Omega = ifelse(adaptive,0.01,1), 
                     delta_Omega = ifelse(adaptive,1e-6,0.01),
                     n_iter = 2000,
                     n_burn_in = 1000, thin_by = 10,
                     progress = T, verbos = T) {
  # some waring messages
  err_no_predictor <- "No predictor supplied.\n\n"
  warr_centering <- "Predictors will be centered.\n\n"

  # check links
  if (!(link %in% c("identity", "probit", "log", "logit"))) {
    stop("Currently only implemented identity (normal), 
      probit (bernoulli) links, log (Poisson) and logit (multi-nomial)")
  }

  # omit NAs
  if (!all(!is.na(data)) & verbos) {
    warning("NA in data will be omitted")
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
    stop("Hyper parameters for B shrinkage must be positive\n\n")
  }
  if (!adaptive) {
    if (is.null(r_beta)) r_beta <- 1
    if (is.null(delta_beta)) delta_beta <- 0.01
    if (verbos & (length(r_beta) > 1 | length(delta_beta) > 1)) {
      warning("Algorithm set to be non-adapive, 
        will take the first entry of hyper prior for B shrinkage\n\n")
    }
    r_beta <- r_beta[1]
    delta_beta <- delta_beta[1]
  }
  else {
    if (is.null(r_beta)) r_beta <- 0.01
    if (is.null(delta_beta)) delta_beta <- 1e-6
    if ((length(r_beta) == 1 & length(delta_beta) == 1)) {
      if (verbos) warning("Algorithm set to be adapive. 
        Assuming all hyper parameters are the same \n\n")
      r_beta <- matrix(r_beta, p, k)
      delta_beta <- matrix(delta_beta, p, k)
    }
    else {
      dims <- c(nrow(r_beta), ncol(r_beta), nrow(delta_beta), ncol(delta_beta))
      prop_dim <- c(p, k, p, k)
      mismatch <- dimname1[dims != prop_dim]
      if (length(mismatch) > 0) {
        errmsg <- paste(
          "Dimension mismatch for hyper prior of B shrinkage: ",
          paste(mismatch, collapse = " "), "\n\n"
        )
        stop(errmsg)
      }
    }
  }
  ## end checking r_beta

  ## checking r_Omega, delta_Omega
  if (is.matrix(r_Omega)) {
    if (verbos) warning("Supplied matrix for hyper parameter for r_Omega, 
      will only take upper triangular part\n\n")
    r_Omega <- c(r_Omega[upper.tri(r_Omega)])
  }

  if (is.matrix(delta_Omega)) {
    if (verbos) warning("Supplied matrix for hyper parameter for delta_Omega, 
      will only take upper triangular part\n\n")
    delta_Omega <- c(delta_Omega[upper.tri(delta_Omega)])
  }

  if (!adaptive) {
    if (is.null(r_Omega)) r_Omega <- 1
    if (is.null(delta_Omega)) delta_Omega <- 0.01
    if (verbos & (length(r_Omega) > 1 | length(delta_Omega) > 1)) {
      warning("Algorithm set to be non-adapive, 
        will take the first entry of hyper prior for Omega shrinkage\n\n")
    }
    r_Omega <- r_Omega[1]
    delta_Omega <- delta_Omega[1]
  }
  else {
    if (is.null(r_Omega)) r_Omega <- 0.01
    if (is.null(delta_Omega)) delta_Omega <- 1e-6
    if ((length(r_Omega) == 1 & length(delta_Omega) == 1)) {
      if (verbos) warning("Algorithm set to be adapive. 
        Assuming all hyper parameters are the same \n\n")
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
  }
  ## end checking Omega hyper parameters

  
  ## Main algorithms

  if (link == "probit") {
    unique_values <- apply(y, 2, function(w) {
      length(unique(w))
    })
    if (!all(unique_values == 2)) {
      stop("Response has multiple unique values, 
        cannot use probit link, do you want logit?\n")
    }

    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    res <- Probit_CAR_LASSO_Cpp(
      y, design, n_iter, n_burn_in,
      thin_by, r_beta, delta_beta,
      r_Omega, delta_Omega, progress
    )
  }

  if (link == "identity") {
    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    if (adaptive) {
      res <- CAR_ALASSO_Cpp(
        y, design,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, progress
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
      res <- Pois_CAR_ALASSO_Cpp(
        y, design,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, progress
      )
    }
    else {
      res <- Pois_CAR_LASSO_Cpp(
        y, design,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, progress
      )
    }
  }

  if (link == "logit") {
    if (verbos) warning("Last response will be used as reference group\n\n")
    if (verbos) cat("Algorithm start...\n\n")
    if (verbos & progress) cat("progress:\n\n")
    if (adaptive) {
      res <- Multinomial_CAR_ALASSO_Cpp(
        y, design,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, progress
      )
    }
    else {
      res <- Multinomial_CAR_LASSO_Cpp(
        y, design,
        n_iter, n_burn_in,
        thin_by, r_beta, delta_beta,
        r_Omega, delta_Omega, progress
      )
    }
  }

  cat("\ndone\n\n")
  res <- lapply(res, coda::mcmc)
  return(res)
}