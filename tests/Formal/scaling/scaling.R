library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(bayesm)

rm(list = ls())

sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")
sourceCpp("./src/CAR-ALASSO.cpp")
source("./R/misc.R")

ks <- c(100,50,25,10,5)
ns <- c( 500,1000 )
ps <- c(10,5)


set.seed(42)

n_exp <- length(ks) * length(ns) * length(ps) * 4
systimet <- (system.time(1+1))
res <- data.frame(matrix(NA,n_exp,9))
colnames(res) <- c("k","n","p","algo",names(systimet))
i_exp <- 1

for(k in ks){
  for(n in ns){
    for(p in ps){
      cat("k=",k,"n=",n,"p=",p,"\n")
      B <- 3*rsparsematrix(k,k,0.2)
      omega <- diag(rgamma(k,3,.1))
      I <- diag(rep(1,k))
      Omega <- t(I-B) %*% omega %*% (I-B)
      Omega <- as.matrix(Omega)
      Sigma <- solve(Omega)
      Design <- 1.0* (matrix(rnorm(n*p,0,1),n,p))
      beta <- as.matrix( rsparsematrix(p,k,0.2))

      mu <-  rnorm(k)

      Xbeta <- Design %*% beta

      Z <- matrix(NA,n,k)

      for( j in 1:n ){
        Z[j,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[j,]+mu),Sigma)
      }
      
      cat("CAR:\n")
      temp <- system.time(CAR_LASSO_Cpp(Z,  Design, n_iter = 1000, 
                          n_burn_in = 100, thin_by = 1, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T))
      temp
      
      res$algo[i_exp] <- "CAR_LASSO"
      res[i_exp, 1:3] <- c(k,n,p)
      res[i_exp, 5:9] <- as.numeric(temp)
      i_exp <- i_exp + 1

      cat("ACAR:\n")
      temp <- system.time(CAR_ALASSO_Cpp(Z,  Design, n_iter = 1000, 
                          n_burn_in = 100, thin_by = 1, 
                          r_beta = 1+0*beta, delta_beta = .01 + 0 * beta,
                          r_Omega = rep(1,.5*(k+1)*k),
                          delta_Omega = rep(.01,.5*(k+1)*k),
                          progress = T))
      temp
      res$algo[i_exp] <- "CAR_ALASSO"
      res[i_exp, 1:3] <- c(k,n,p)
      res[i_exp, 5:9] <- as.numeric(temp)
      i_exp <- i_exp + 1

      cat("Glasso:\n")
      temp <- system.time(Graphical_LASSO_Cpp(Z, 1000, 
                          100, 10, 1, .01, T))

      temp
      res$algo[i_exp] <- "GLASSO"
      res[i_exp, 1:3] <- c(k,n,p)
      res[i_exp, 5:9] <- as.numeric(temp)
      i_exp <- i_exp + 1

      cat("multireg:\n")
      temp <- system.time(sample_multireg <- lapply(1:1100,function(i) 
                rmultireg(Z,cbind(1,Design),
                0*rbind(1,beta),diag(0.001,p+1,p+1),3,3*diag(2,k,k))))

      temp
      res$algo[i_exp] <- "multireg"
      res[i_exp, 1:3] <- c(k,n,p)
      res[i_exp, 5:9] <- as.numeric(temp)
      i_exp <- i_exp + 1
      
      write.csv(res,"./tests/Formal/scaling/scalingres.csv")
    }
  }
}