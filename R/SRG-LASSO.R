SRGlasso <- function(formula,# a double sided formula needed, e.g. x+y~a+b
                     data,link = "identity",
                     n_iter = 2000, 
                     n_burn_in = 1000, thin_by = 10, 
                     r_beta = 1, delta_beta = .1,
                     r_Omega = 1,delta_Omega = .1,
                     progress = T,verbos = T){
  require(Matrix)
  require(Rcpp)
  require(RcppArmadillo)
  require(RcppProgress)
  require(coda)
  
  sourceCpp("./src/Probit-SRG-LASSO.cpp")
  sourceCpp("./src/SRG-LASSO.cpp")
  sourceCpp("./src/graphical-LASSO.cpp")
  sourceCpp("./src/Probit-Graphical-LASSO.cpp")
  
  
  if(!(link %in% c("identity","probit"))){
    stop("Currently only implemented identity (normal) and probit (bernoulli) links")
  }
  
  if(!all(!is.na(data))){
    warning("NA in data will be omitted")
    data <- na.omit(data)
  }
  
  design <- model.matrix(formula,data)
  design <- design[,colnames(design)!="(Intercept)"] # design matrix, no intercept
  design <- as.matrix(design)
  
  
  response <- formula
  response[[3]] <- formula[[2]]
  y <- model.matrix(response,data) # response matrix.
  y <- y[,colnames(y)!="(Intercept)"]
  rm(response)
  
  if(link == "probit"){
    unique_values <- apply(y,2,function(w){length(unique(w))})
    if(!all(unique_values==2)){
      stop("Response has multiple unique values, cannot use probit link\n")
    }
    
    if(ncol(design)==0){
      if(verbos) cat("No predictor other than intercept\n\nAlgorithm start...\n\n")
      cat("progress:\n\n")
      Res <- Probit_Graphical_LASSO_Cpp(y,n_iter,n_burn_in,
                                        thin_by,r_Omega,delta_Omega,progress)
    }
    else{
      if(verbos){
        cat("Predictors will be centered and intercept will be included...\n\n")
      }
      design <- apply(design,2,function(w){(w-mean(w))/sd(w)})
      if(verbos) cat("Algorithm start...\n\n")
      cat("progress:\n\n")
      Res <- Probit_SRG_LASSO_Cpp(y,design,n_iter,n_burn_in,
                                  thin_by,r_beta,delta_beta,
                                  r_Omega,delta_Omega,progress)
    }
    
    
  }
  if(link == "identity"){
    
    if(ncol(design)==0){
      if(verbos) cat("No predictor other than intercept, will use graphical-LASSO\n\nAlgorithm start...\n\n")
      cat("progress:\n\n")
      Res <- Graphical_LASSO_Cpp(y,n_iter,n_burn_in,
                                        thin_by,r_Omega,delta_Omega,progress)
    }
    else{
      if(verbos){
        cat("Predictors will be centered and intercept will be included...\n\n")
      }
      
      design <- apply(design,2,function(w){(w-mean(w))/sd(w)})
      if(verbos) cat("Algorithm start...\n\n")
      cat("progress:\n\n")
      Res <- SRG_LASSO_Cpp(y,design,n_iter,n_burn_in,
                                  thin_by,r_beta,delta_beta,
                                  r_Omega,delta_Omega,progress)
    }
    
    
  }
  
  cat("\ndone\n\n")
  Res <- lapply(Res,coda::mcmc)
  return(Res)
}