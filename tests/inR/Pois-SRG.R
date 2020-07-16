library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 11
n = 8000
p = 1

sourceCpp("./src/Pois-SRG-LASSO.cpp")


set.seed(43)
B <- 0.3*rsparsematrix(k,k,0.3)
omega <- diag(rgamma(k,15,.5))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
#diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)
image(Omega)

Sigma <- solve(Omega)

Design <- matrix(rnorm(n*p,0,1),n,p)
Design <- (Design-mean(Design))/sd(Design)
colnames(Design) <- paste0("x",1:p)


beta <- matrix(rnorm(p*k,0,.5),p,k)
#beta[sample(p*k,floor(0.3*p*k))] = 0

mu <- rnorm(k,2,.5)
#mu

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)
Y <- Z

for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Xbeta[i,]+mu,Sigma)
  Y[i,] <- rpois(k,exp(Z[i,]))
}

Z_back <- 1*Z
Z <- 1*Z_back

sum(is.na(Z))
update_Z_helper_Pois_reg(Z,Y,Design,mu,beta,Omega,k,p,n,700,3,64)
sum(is.na(Z))

kkk <- Pois_SRG_LASSO_Cpp(Y, 
                   Design, 
                   5000, # how many iterations?
                   200, # burn in
                   1, # thinning?
                   1, # prior on lambda of beta
                   .01,
                   1,
                   .01,
                   100, 5, 64, # ars parameters
                   T)



f <- function(x,y,mu,sigma2){
  y*x-exp(x)-0.5*(x-mu)^2/sigma2
}

fprime <- function(x,y,mu,sigma2){
  y-exp(x)-(x-mu)/sigma2
  
}

y=1
m = 4
ars::ars(1,f,fprime,y=y,mu=-5.35,sigma2 = 0.5975,x = (log(y+.1)-5.35)/2 + ((1:m-1)-(m/2))*(4*(y+abs(-5.35)))/m,m=m)



