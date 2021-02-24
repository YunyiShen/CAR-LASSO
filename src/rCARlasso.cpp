// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "CAR_LASSO_helper.h"

//' Block Gibbs sampler for CAR-LASSO
//' @description \strong{This function is for advanced users to build their own sampler use CARlasso as core.} It will take one round of Gibbs sampler of CAR-LASSO model. Be aware that the function is a `void` function implemented in C++, and all updated parameters e.g. Omega will be manipulate directly in memory to save space. Users should manage to do their own work to save the state. Also be aware that R uses shallow copy by default, which means one cannot save the state by simply give it to another object e.g. first `Omega_old <- Omega_curr` then update `Omega_curr`, `Omega_old` will also change. \strong{This function will NOT check dimensions of input.} Below we assule n data, k responses and p predictors.
//' @param Z_curr the current (latent) normal data, should be n*k. Will not be changed
//' @param design the design matrix, should be n*p. Will not be changed
//' @param lambda2_beta the current shrinkage parameter of regression coefficients, should be a scalar of type `double`. Will be updated
//' @param tau2_curr the current latent scale parameter in the normal mixture representation of Laplace, for regression coefficients, should be a vector with p*k entries. Will be updated.
//' @param beta_curr the current regression coefficients, should be a matrix sized p*k (p row and k columns). Will be updated.
//' @param lambda_Omega the current shrinkage parameter for Omega, should be a scalar of tyoe `double`. Will be updated.
//' @param Omega_curr the current Omega matrix, should be a matrix of size k*k. Will be updated.
//' @param mu_curr the current mu, intercept, should be a vector of size k. Will be updated.
//' @param r_beta hyperprior's parameter of shrinkage for regression coefficients, should be a scalar of type 'double' and positive. Will not be updated.
//' @param delta_beta hyperprior's parameter of shrinkage for regression coefficients, should be a scalar of type 'double' and positive. Will not be updated.
//' @param r_Omega hyperprior's parameter of shrinkage for precision Omega, should be a scalar of type 'double' and positive. Will not be updated.
//' @param delta_beta hyperprior's parameter of shrinkage for rprecision Omega, should be a scalar of type 'double' and positive. Will not be updated.
//' @param k integer, number of responses
//' @param p integer, number of predictors
//' @param n integer, number of data points
//' @return Again this is a `void` function and will not return anything. All update happened in memory dirrectly. 


// [[Rcpp::export]]
void rCARlasso(const arma::mat &Z_curr,
               const arma::mat &design, //design matrix
               double &lambda2_beta,    // lambda2, the shrinkage for beta, double
               arma::vec &tau2_curr,    // tau2 for beta, should be (k*p) entries vector
               arma::mat &beta_curr,    // current value of beta, should be p*k matrix
               double &lambda_Omega,    // lambda, the shrinkage for Omega, double
               arma::mat &Omega_curr,   // current value of Omega
               arma::vec &mu_curr,      // current value of mu
               const double r_beta,     // prior on lambda of beta
               const double delta_beta,
               const double r_Omega, // prior on lambda of Omega
               const double delta_Omega,
               int k, int p, int n)
{

    double Omega_r_post = (r_Omega + (k * (k + 1) / 2));
    double Omega_delta_post;
    Omega_delta_post = (delta_Omega + sum(sum(abs(Omega_curr))) / 2);
    lambda_Omega = randg<double>(distr_param(Omega_r_post, 1 / Omega_delta_post));
    //Update betas:
    beta_curr = update_car_beta_helper(Z_curr, design, mu_curr,
                                       tau2_curr, Omega_curr,
                                       k, p, n);

    // update Omega

    update_car_Omega_helper(Omega_curr, Z_curr, design,
                            mu_curr, beta_curr,
                            lambda_Omega,
                            k, p, n);

    // Update mu

    mu_curr = update_car_mu_helper(Z_curr, design, beta_curr,
                                   Omega_curr,
                                   k, p, n);

    // Update tau2 for beta
    tau2_curr = update_car_tau2_helper(beta_curr, lambda2_beta,
                                       Omega_curr, k, p, n);

    // Update lambda_beta

    lambda2_beta = randg<double>(distr_param(r_beta + k * p, 1 / (delta_beta + sum(tau2_curr) / 2)));

    return;
}