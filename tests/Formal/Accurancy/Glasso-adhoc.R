library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)
library(bayesm)



sourceCpp("./src/graphical-LASSO.cpp")
source("./tests/Formal/Accurancy/misc.R")
source("./tests/Formal/Accurancy/Graph_generator.R")
settings_link <- "./tests/Formal/Accurancy/settings/"

ks <- c(30,100)
ss <- c(.2,.5)
ps <- c(5,10)
nrep <- 50

thr <- c( 1e-2 ,5e-3, 1e-3 )

for(k in ks){
    n <- 50 * (k==30) + 200 * (k==100)
    for(p in ps){
        for(s in ss){

            beta <- read.csv(paste0(
                settings_link, "beta_k", k, "_p", p , "_s",s, ".csv"
            ))
            beta <- as.matrix(beta)
            design <- read.csv(paste0(
                settings_link, "design_p",p,"_n",n,".csv"
            ))
            design <- as.matrix(design)
            
        }
    }

}
