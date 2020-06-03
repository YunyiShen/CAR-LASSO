// [[Rcpp::depends(RcppArmadillo)]]
#include <tgmath.h>
#include <RcppArmadillo.h> 
#include <iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"

// [[Rcpp::export]]
arma::vec update_tau2_helper(const arma::mat & beta,
                             const double & lambda2,
                             const arma::mat & Omega,
                             int k, int p, int n){
  arma::vec betavec = vectorise(beta);
  arma::vec invtau2(k*p);
  
  double detSigma = 1/det(Omega);
  detSigma = pow(detSigma,1/(2*k));
  
  
  
  arma::vec mu_prime = sqrt(lambda2*detSigma/(betavec%betavec));
  for(int i = 0 ; i < k*p ; ++i){
    invtau2(i) =  rinvGau(mu_prime(i),lambda2);
  }
  return(1/invtau2);
}

