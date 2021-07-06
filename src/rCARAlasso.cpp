// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "CAR_LASSO_helper.h"
#include "CAR_ALASSO_helper.h"

//' @title Block Gibbs sampler for adaptive CAR-LASSO
//'
//' @description \strong{This function is for advanced users to build their own sampler use adaptive CARlasso as core.} It will execute one round of Gibbs sampler of adaptive CAR-LASSO model. Be aware that the function is a `void` function implemented in C++, and all updated parameters e.g. Omega will be manipulate directly in memory to save space. Users should manage to do their own work to save the state. Also be aware that R uses shallow copy by default, which means one cannot save the state by simply give it to another object e.g. first `Omega_old <- Omega_curr` then update `Omega_curr`, `Omega_old` will also change. \strong{This function will NOT check dimensions of input.} Below we assume n samples, k responses and p predictors.
//' @param Z_curr the current (latent) normal Z_curr, should be n*k. Will not be changed
//' @param design the design matrix, should be n*p. Will not be changed
//' @param lambda2_beta the current shrinkage parameter of regression coefficients, should be a vector with p*k entries. Will be updated
//' @param tau2_curr the current latent scale parameter in the normal mixture representation of Laplace, for regression coefficients, should be a vector with p*k entries. Will be updated.
//' @param beta_curr the current regression coefficients, should be a matrix sized p*k (p row and k columns). Will be updated.
//' @param lambda_Omega the current shrinkage parameter for Omega, should be a vector with k*(k-1)/2 entries. Will be updated.
//' @param Omega_curr the current Omega matrix, should be a matrix of size k*k. Will be updated.
//' @param mu_curr the current mu, intercept, should be a vector of size k. Will be updated.
//' @param r_beta hyperprior's parameter of shrinkage for regression coefficients, should be a scalar of type 'double' and positive. Will not be updated.
//' @param delta_beta hyperprior's parameter of shrinkage for regression coefficients, should be a scalar of type 'double' and positive. Will not be updated.
//' @param r_Omega hyperprior's parameter of shrinkage for precision Omega, should be a scalar of type 'double' and positive. Will not be updated.
//' @param delta_Omega hyperprior's parameter of shrinkage for rprecision Omega, should be a scalar of type 'double' and positive. Will not be updated.
//' @param lambda_diag shrinkage parameter of the diagonal of Omega, should be a vector of size k, should be non-negative. Will not be updated
//' @param k integer, number of responses
//' @param p integer, number of predictors
//' @param n integer, number of Z_curr points
//' @return Again this is a `void` function and will not return anything. All update happened in memory directly. 
// [[Rcpp::export]]
void rCARAlasso_(const arma::mat &Z_curr,
               const arma::mat &design, //design matrix
               arma::vec &lambda2_beta,    // lambda2, the shrinkage for beta, vec size k*p
               arma::vec &tau2_curr,    // tau2 for beta, should be (k*p) entries vector
               arma::mat &beta_curr,    // current value of beta, should be p*k matrix
               arma::vec &lambda_Omega,    // lambda, the shrinkage for Omega, double
               arma::mat &Omega_curr,   // current value of Omega
               arma::vec &mu_curr,      // current value of mu
               const arma::mat & r_beta,     // prior on lambda of beta
               const arma::mat & delta_beta,
               const arma::vec & r_Omega, // prior on lambda of Omega
               const arma::vec & delta_Omega,
               const arma::vec & lambda_diag,
               int k, int p, int n)
{
    // Update lambda_Omega
    update_car_lambda_Omega_adp_helper(lambda_Omega,
                                       Omega_curr,
                                       r_Omega,
                                       delta_Omega);
    
    //Update betas:
    beta_curr = update_car_beta_helper(Z_curr, design, mu_curr,
                                   tau2_curr, Omega_curr, 
                                   k, p, n);
    
    // flaged
    // update Omega
    
    update_car_Omega_adp_helper(Omega_curr, Z_curr, design, 
                                     mu_curr, beta_curr,
                                     lambda_Omega,
                                     lambda_diag,
                                     k, p, n);
    
    //Rcout << "flag" <<endl;
    
    // Update mu
    
    mu_curr = update_car_mu_helper(Z_curr,design,beta_curr,
                               Omega_curr, 
                               k, p, n);
    
    
  
    // Update tau2 for beta
    tau2_curr = update_car_tau2_adp_helper(beta_curr,lambda2_beta,
                                   Omega_curr,k,p,n);
    
    
    // Update lambda_beta
    
    update_car_lambda2_beta_adp_helper(lambda2_beta,tau2_curr,
                                        r_beta,delta_beta,k, p);


    return; 


}

