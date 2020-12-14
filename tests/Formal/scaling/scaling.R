library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(bayesm)
library(ggplot2)

rm(list = ls())

sourceCpp("./src/CAR-LASSO.cpp")
sourceCpp("./src/SRG-LASSO.cpp")
sourceCpp("./src/graphical-LASSO.cpp")
sourceCpp("./src/CAR-ALASSO.cpp")
source("./R/misc.R")

ks <- c(100,50,25,10,5)
ns <- c( 500,1000 )
ps <- c(10,5)


set.seed(42)

n_exp <- length(ks) * length(ns) * length(ps) * 5
systimet <- (system.time(1+1))
res <- data.frame(matrix(NA,n_exp,9))
colnames(res) <- c("k","n","p","algo",names(systimet))
i_exp <- 1

for(k in ks){
  for(p in ps){
    
      
    B <- 3*rsparsematrix(k,k,0.2)
    omega <- diag(rgamma(k,3,.1))
    I <- diag(rep(1,k))
    Omega <- t(I-B) %*% omega %*% (I-B)
    Omega <- as.matrix(Omega)
    Sigma <- solve(Omega)
      
    beta <- as.matrix( rsparsematrix(p,k,0.2))

    mu <-  rnorm(k)
      
    for(n in ns){  
      cat("k=",k,"n=",n,"p=",p,"\n")
      Design <- 1.0* (matrix(rnorm(n*p,0,1),n,p))
      Xbeta <- Design %*% beta
      Z <- matrix(NA,n,k)

      for( j in 1:n ){
        Z[j,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[j,]+mu),Sigma)
      }
      
      cat("SRG:\n")
      temp <- system.time(SRG_test <- SRG_LASSO_Cpp(Z,  Design, n_iter = 1000, 
                          n_burn_in = 100, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T))
      
      temp
      
      res$algo[i_exp] <- "SRG_LASSO"
      res[i_exp, 1:3] <- c(k,n,p)
      res[i_exp, 5:9] <- as.numeric(temp)
      i_exp <- i_exp + 1
      write.csv(res,"./tests/Formal/scaling/scalingres.csv")

      cat("CAR:\n")
      temp <- system.time(CAR_LASSO_Cpp(Z,  Design, n_iter = 1000, 
                          n_burn_in = 100, thin_by = 10, 
                          r_beta = 1, delta_beta = .01,
                          r_Omega = 1,delta_Omega = .01,
                          progress = T))
      temp
      
      res$algo[i_exp] <- "CAR_LASSO"
      res[i_exp, 1:3] <- c(k,n,p)
      res[i_exp, 5:9] <- as.numeric(temp)
      i_exp <- i_exp + 1
      write.csv(res,"./tests/Formal/scaling/scalingres.csv")

      cat("ACAR:\n")
      temp <- system.time(CAR_ALASSO_Cpp(Z,  Design, n_iter = 1000, 
                          n_burn_in = 100, thin_by = 10, 
                          r_beta = 1+0*beta, delta_beta = .01 + 0 * beta,
                          r_Omega = rep(1,.5*(k+1)*k),
                          delta_Omega = rep(.01,.5*(k+1)*k),
                          progress = T))
      temp
      res$algo[i_exp] <- "CAR_ALASSO"
      res[i_exp, 1:3] <- c(k,n,p)
      res[i_exp, 5:9] <- as.numeric(temp)
      i_exp <- i_exp + 1
      write.csv(res,"./tests/Formal/scaling/scalingres.csv")

      cat("Glasso:\n")
      temp <- system.time(Graphical_LASSO_Cpp(Z, 1000, 
                          100, 10, 1, .01, T))

      temp
      res$algo[i_exp] <- "GLASSO"
      res[i_exp, 1:3] <- c(k,n,p)
      res[i_exp, 5:9] <- as.numeric(temp)
      i_exp <- i_exp + 1
      write.csv(res,"./tests/Formal/scaling/scalingres.csv")
      
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


res = read.csv("scalingres.csv", header=TRUE)

res$n <- paste0("sample_size:",res$n)
#res$p <- paste0("predictors:", res$p)
res$p <- paste0(res$p, " predictors")
res <- within(res, p<-factor(p, levels=c("5 predictors", "10 predictors")))
res <- within(res, n<-factor(n, levels=c("sample_size:500", "sample_size:1000")))
res <- within(res, algo<-factor(algo))
res <- within(res, algo<-factor(algo, levels=c("CAR_LASSO", "CAR_ALASSO", "SRG_LASSO", "GLASSO", "multireg")))
res <- within(res, algo<-factor(algo, labels=c("CAR-LASSO", "CAR-ALASSO", "SRG-LASSO", "GLASSO", "multireg")))



cbbPalette <- c(
                "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7")

ggplot(res,
       aes(k, elapsed, col=algo)) + 
  scale_colour_manual(values=cbbPalette) +
  geom_point() +
  geom_line()+
  xlab("number of nodes")+
  ylab("CPU time(s)/1.1k samples")+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  labs(color = " ") + theme_bw() + 
  theme(legend.position="top", legend.box = "vertical") + 
  scale_fill_brewer()+
  theme(text = element_text(size=14), 
        legend.text=element_text(size=8),
        #axis.text.x = element_text(angle=0,size = 12),
        plot.margin = margin(.15, .15, .15, .15, "cm"))+
#  facet_grid(p~n,labeller = label_parsed, scales = 'fixed')
  facet_grid(p~n,scales = 'fixed')

ggsave("scaling_test2.pdf",width = 6,height = 6, unit = "in")
#ggsave("scaling_test.pdf",width = 6,height = 5, unit = "in")
#ggsave("scaling_test.png",width = 6,height = 5, unit = "in", dpi = 500)
