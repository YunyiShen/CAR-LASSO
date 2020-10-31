// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"
#include "CAR_LASSO_helper.h"
#include "CAR_LASSO_randeff_helper.h"


// helper function for updating random effect matrix, so that the conditional mean was design_r * res

// [[Rcpp::export]]
arma::mat update_car_nu_helper(const arma::mat & data,
                               const arma::mat & design, // design matrix
                               const arma::mat & design_r, // design mat for random effect
                               const arma::mat & beta,
                               const arma::vec & mu, // grand mean
                               const arma::vec & xi, // the precision vector of random effect, should be size of k
                               const arma::mat & Omega,
                               int k, int pr, int n){

  arma::mat k_1s(pr,1,fill::zeros);
  k_1s += 1;
  arma::mat Q_nu(k*pr,k*pr,fill::zeros);
  
  arma::mat XtX = design_r.t() * design_r ;

  arma::mat Sigma = inv_sympd(Omega);

  arma::mat mu_nu_mat = data*Omega - design * beta;
  mu_nu_mat.each_row() -= mu.t();

  arma::vec mu_nu = vectorise(mu_nu_mat);
  
  Q_nu = kron(Sigma,XtX); // precision matrix of beta, more intuitive way was sum_i X_i^TSigmaX_i, but kron is faster
  arma::mat perc_random = kron(xi,k_1s);
  Q_nu.diag() += perc_random;


  arma::mat res(size(mu_nu),fill::randn);

  arma::mat chol_Q = arma::chol(Q_nu);
  res = arma::solve(chol_Q,res) + arma::solve(Q_nu,mu_nu);

  res = reshape(res,pr,k);
  
  return(res);
  
}

// [[Rcpp::export]]
void update_xi_helper(arma::vec xi,
                      const arma::mat nu,
                      const double & alpha,
                      const double & beta, 
                      int k, int pr){
    double alpha_post = alpha + pr/2;
    double beta_post;
    for(int i = 0; i < k; ++i){
        beta_post = beta + arma::as_scalar( arma::sum(nu.col(i) * nu.col(i),0) )/2;
        xi(i) = randg<double>( distr_param(alpha_post,1/beta_post) );
    }
    return;
}