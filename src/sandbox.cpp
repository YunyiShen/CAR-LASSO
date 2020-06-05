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
#include "SRG_LASSO_helper.h"
#include "Truncated_Normal_helper.h"

// [[Rcpp::export]]
void update_Z_helper(arma::mat & Z_curr, // persumably large, thus will not copy
                     const arma::mat & data, 
                     const arma::mat & design,
                     const arma::vec & mu_curr,
                     const arma::mat & beta_curr,
                     const arma::mat & Omega_curr,
                     int k, int p, int n){
  arma::mat mu_Zmat = design * beta_curr;
  mu_Zmat.each_row() += mu_curr.t(); // calculate the expectation of latent
  arma::mat Sigma_Z = Omega_curr;
  Sigma_Z.diag() += 1;
  Sigma_Z = inv_sympd(Sigma_Z);
  
  arma::vec mu_Zi;
  
  arma::vec y_star(k,fill::zeros); // latent y (y|Z~ N(Z,1)), bascally the probit trnasformation
  
  for(int i = 0 ; i < n ; ++i){
    for(int j = 0 ; j < k ; ++j){
      y_star(i) = rtn1(Z_curr(i,j),1,
             data(i,j) == 1 ? 0 : -INFINITY, // if data(i,j)=1, then y_star >= 0
             data(i,j) == 0 ? 0 : INFINITY); // if data(i,j)=0 y_star<0
    }
    mu_Zi = Sigma_Z * (Omega_curr * trans(mu_Zmat.row(i))+y_star);
    Z_curr.row(i) = trans( mvnrnd(mu_Zi,Sigma_Z));
  }
}

