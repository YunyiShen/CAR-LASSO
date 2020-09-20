// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
#include <progress.hpp>
#include <progress_bar.hpp>

#include "Truncated_Normal_helper.h"
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
      y_star(j) = rtn1(Z_curr(i,j),1,
             data(i,j) == 1 ? 0 : -INFINITY, // if data(i,j)=1, then y_star >= 0
             data(i,j) == 0 ? 0 : INFINITY); // if data(i,j)=0 y_star<0
    }
    mu_Zi = Sigma_Z * (Omega_curr * trans(mu_Zmat.row(i))+y_star);
    Z_curr.row(i) = trans( mvnrnd(mu_Zi,Sigma_Z));
  }
}


void update_Z_helper_CAR(arma::mat & Z_curr, // persumably large, thus will not copy
                     const arma::mat & data, 
                     const arma::mat & design,
                     const arma::vec & mu_curr,
                     const arma::mat & beta_curr,
                     const arma::mat & Omega_curr,
                     int k, int p, int n){
  arma::mat mu_Zmat = design * beta_curr;
  mu_Zmat.each_row() += mu_curr.t(); // calculate the expectation of latent

  //mu_Zmat = arma::trans(arma::solve(Omega_curr,mu_Zmat.t()));// the coupled structure

  arma::mat Sigma_Z = Omega_curr;
  Sigma_Z.diag() += 1;
  Sigma_Z = inv_sympd(Sigma_Z);
  
  arma::vec mu_Zi;
  
  arma::vec y_star(k,fill::zeros); // latent y (y|Z~ N(Z,1)), bascally the probit trnasformation
  
  for(int i = 0 ; i < n ; ++i){
    for(int j = 0 ; j < k ; ++j){
      y_star(j) = rtn1(Z_curr(i,j),1,
             data(i,j) == 1 ? 0 : -INFINITY, // if data(i,j)=1, then y_star >= 0
             data(i,j) == 0 ? 0 : INFINITY); // if data(i,j)=0 y_star<0
    }
    mu_Zi = Sigma_Z * (trans(mu_Zmat.row(i))+y_star);
    Z_curr.row(i) = trans( mvnrnd(mu_Zi,Sigma_Z));
  }
}



void update_Z_graphical_helper(arma::mat & Z_curr, // persumably large, thus will not copy
                               const arma::mat & data, 
                               const arma::vec & mu_curr,
                               const arma::mat & Omega_curr,
                               int k, int n){
  arma::mat Sigma_Z = Omega_curr;
  Sigma_Z.diag() += 1;
  Sigma_Z = inv_sympd(Sigma_Z);
  
  arma::vec mu_Zi;
  
  arma::vec y_star(k,fill::zeros); // latent y (y|Z~ N(Z,1)), bascally the probit trnasformation
  
  for(int i = 0 ; i < n ; ++i){
    for(int j = 0 ; j < k ; ++j){
      y_star(j) = rtn1(Z_curr(i,j),1,
             data(i,j) == 1 ? 0 : -INFINITY, // if data(i,j)=1, then y_star >= 0
             data(i,j) == 0 ? 0 : INFINITY); // if data(i,j)=0 y_star<0
    }
    mu_Zi = Sigma_Z * (Omega_curr * mu_curr +y_star);
    Z_curr.row(i) = trans( mvnrnd(mu_Zi,Sigma_Z));
  }
}
