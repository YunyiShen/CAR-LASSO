library(Matrix)
library(RcppArmadillo)
library(Rcpp)
library(RcppProgress)

rm(list = ls())

k = 11
n = 8000
p = 2


sourceCpp("./src/CAR_LASSO_helper.cpp")

B <- rsparsematrix(k,k,0.3)
omega <- diag(rgamma(k,10,.3))
I <- diag(rep(1,k))
Omega <- t(I-B) %*% omega %*% (I-B)
#diag(Omega) <- diag(Omega) + k
Omega <- as.matrix(Omega)
image(Omega)

Sigma <- solve(Omega)

Design <- matrix(rnorm(n*p,0,1),n,p)
Design <- (Design-mean(Design))/sd(Design)
colnames(Design) <- paste0("x",1:p)


beta <- matrix(rnorm(p*k,5,1),p,k)
#beta[sample(p*k,floor(0.3*p*k))] = 0

mu <-  rnorm(k)
#mu

Xbeta <- Design %*% beta

Z <- matrix(NA,n,k)


for( i in 1:n ){
  Z[i,] <- MASS::mvrnorm(1,Sigma %*% (Xbeta[i,]+mu),Sigma)
  #Y[i,] <- 1 * ((rnorm(k,Z[i,],1))>0)
}


n_iter = 500
n_burnin = 100
## beta

beta_save = list()
beta_temp = beta
lambda_temp = 0.2
for(i in 1:(n_iter+n_burnin)){
  tau2 = update_car_tau2_helper(beta_temp,lambda_temp,Omega,k,p,n)
  beta_temp = update_car_beta_helper(Z,Design,mu,tau2,Omega,k,p,n)
  lambda_temp = rgamma(1,1+k*p,sum(tau2)/2+.1)
  cat(lambda_temp,"\n")
  if(i>n_burnin){
    beta_save[[i-n_burnin]] = beta_temp
  }
  #cat(i,"\n")
}
Reduce("+",beta_save)/n_iter

## mu 

mu_save = list()

for(i in 1:(n_iter+n_burnin)){
  mu_temp = update_car_mu_helper(Z,Design,beta,Omega,k,p,n)
  
  if(i>n_burnin){
    mu_save[[i-n_burnin]] = mu_temp
  }
  cat(i,"\n")
}
c(Reduce("+",mu_save)/n_iter)



Omega_save = list()

Omega_going_to_be_saved <- solve(cov(Z))

for(i in 1:(n_iter+n_burnin)){
  
  
  if(i>n_burnin & (i-n_burnin)%%1000 == 0){
    temp = Omega_going_to_be_saved
    temp[1,1] = 0
    temp[1,1] = Omega_going_to_be_saved[1,1]
    Omega_save[[(i-n_burnin)/1000]] = temp
  }
  update_car_Omega_helper(Omega_going_to_be_saved,Z,
                          Design,mu,beta,.1,k,p,n)
  cat(i,"\n")
}
mean_Omega = Reduce("+",Omega_save)/(n_iter/1000)
image(mean_Omega)
image(Omega)

hist((mean_Omega-Omega)/Omega)
hist((solve(cov(Z)))/Omega)

require(coda)
acf(mcmc(sapply(Omega_save,"[",1)))


## combine
Omega_save = list()
beta_save = list()
mu_save = list()

Omega_going_to_be_saved <- solve(cov(Z))

for(i in 1:(n_iter+n_burnin)){
  Omega_curr=Omega_going_to_be_saved
  update_car_Omega_helper(Omega_going_to_be_saved,Z,
                          Design,mu,beta,8,k,p,n)
  
  
  
  if(i>n_burnin){
    temp = Omega_going_to_be_saved
    temp[1,1] = 0
    temp[1,1] = Omega_going_to_be_saved[1,1]
    Omega_save[[i-n_burnin]] = temp
  }
  cat(i,"\n")
}




