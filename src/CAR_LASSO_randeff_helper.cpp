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
                               const arma::mat & membership, // membership matrix for random precision, this should be something similar to a design matrix, that tells the algorithm which precision correspond to which random effect, important for multiple memberships, then the precision will be membership * xi
                               const arma::mat & beta,
                               const arma::vec & mu, // grand mean
                               const arma::mat & xi, // the precision matrix of random effect, should be m rows k column, m is how many different distributions we have
                               const arma::mat & Omega,
                               int k, int pr, int n){

  arma::mat Q_nu(k*pr,k*pr,fill::zeros);
  
  arma::mat XtX = design_r.t() * design_r ;

  arma::mat Sigma = inv_sympd(Omega);

  arma::mat mu_nu_mat = data*Omega - design * beta;
  mu_nu_mat.each_row() -= mu.t();
  mu_nu_mat = design_r.t() * (mu_nu_mat * Sigma);

  arma::vec mu_nu = vectorise(mu_nu_mat);
  
  Q_nu = kron(Sigma,XtX); // precision matrix of beta, more intuitive way was sum_i X_i^TSigmaX_i, but kron is faster
  arma::mat perc_random = membership * xi;
  //Rcout << perc_random << endl;
  Q_nu.diag() += vectorise(perc_random);
  //Rcout << Q_nu << endl;
  //Rcout << mu_nu << endl;

  arma::mat res(size(mu_nu),fill::randn);

  arma::mat chol_Q = arma::chol(Q_nu);
  res = arma::solve(chol_Q,res) + arma::solve(Q_nu,mu_nu);

  res = reshape(res,pr,k);
  
  return(res);
  
}

// [[Rcpp::export]]
void update_xi_helper(arma::mat & xi,
                      const arma::mat & nu,
                      const arma::mat & membership,
                      const double & alpha,
                      const double & beta, 
                      int k, int pr, int m){
    double alpha_post; 
    double beta_post;
    int nn; // number of latent variable for that member
    arma::vec latent_m; //vector of latent variable of mth membership
    for(int i = 0; i < m; ++i){
        for (int j = 0; j < k; ++j)
        {
            nn = as_scalar(sum(membership.col(i))); // number of latent var for node i, member j 
            alpha_post = alpha + nn/2;
            latent_m = membership.col(i) % nu.col(j);// this vector was 0 for non members and nu for members
            beta_post = beta + arma::as_scalar( arma::sum(  latent_m % latent_m) )/2;
            xi(i,j) = randg<double>( distr_param(alpha_post,1/beta_post) );
        }
    }
    return;
}

// [[Rcpp::export]]
void get_data_centered(arma::mat & centered_data,
                       const arma::mat & data,
                       const arma::mat & design_r,
                       const arma::mat & nu,
                       const arma::mat & Omega){
    arma::mat Sigma = inv(Omega);
    centered_data = data - (design_r * nu) * Sigma;
    return;
}