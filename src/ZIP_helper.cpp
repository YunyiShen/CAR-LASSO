// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
#include "Truncated_Normal_helper.h"

/*
 * Calulate latent variables of the Pois part
 */
void ZIP_update_Pois_latent_helper(mat & Pois_latent,
                                   const mat & data,
                                   const mat & Z_pois,
                                   const mat & Z_Prob,
                                   int k, int p, int n){
  for(int i = 0; i < n ; ++i){
    for(int j = 0 ; j < k ; ++j){
      Pois_latent(i,j) = Z_prob(i,j)<=0 ? R::rpois(exp(Z_pois(i,j))) : data(i,j);
    }
  }
  return;
}


void ZIP_update_Z_Prob_helper(mat & Z_prob,
                              const mat & data,
                              const arma::mat & design,
                              const arma::vec & mu_curr,
                              const arma::mat & beta_curr,
                              const arma::mat & Omega_curr,
                              const mat & Pois_latent,
                              int k, int p, int n){
  arma::mat mu_Zmat = design * beta_curr;
  mu_Zmat.each_row() += mu_curr.t(); // calculate the expectation of latent
  arma::mat Sigma_Z = Omega_curr;
  Sigma_Z.diag() += 1;
  Sigma_Z = inv_sympd(Sigma_Z);
  
  arma::vec mu_Zi;
  arma::vec y_star(k,fill::zeros); // latent y (y|Z~ N(Z,1)), bascally the probit trnasformation
  
  int temp = 0; // 0-inflate or not
  
  
  for(int i = 0; i < n ; ++i){
    for(int j = 0 ; j < k ; ++j){
      if(data(i,j)==0){
        temp = Pois_latent(i,j)==0 ? (1.0 * R::rnorm(Z_prob(i,j),1) <= 0) : 0 ;
      }
      else{
        temp = 1 ;
      }
      y_star(j) = rtn1(Z_prob(i,j),1,
             temp == 1 ? 0 : -INFINITY, // if data(i,j)=1, then y_star >= 0
             temp == 0 ? 0 : INFINITY); // if data(i,j)=0 y_star<0
    }
    mu_Zi = Sigma_Z * (Omega_curr * trans(mu_Zmat.row(i))+y_star);
    Z_prob.row(i) = trans( mvnrnd(mu_Zi,Sigma_Z));
  }
  return;
}

