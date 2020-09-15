#ifndef ARS_POIS_HELPER_H
#define ARS_POIS_HELPER_H
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
//#include "Error.h"
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <iostream>
#include <climits>
#include <cmath>

void sample_(int *iwv,      double *rwv, 
             int *ifault, int l, int w,// which node
             arma::mat & Z_curr,
             const arma::mat & mu_Z,
             const arma::mat & Sigma_Z,// this is Sigma (cov) not Omega (percision)
             const arma::mat & y,
             int k, int p, int n_sample);


void update_Z_helper_Pois(arma::mat & Z_curr,
                          const arma::mat & mu_Z,
                          const arma::mat & Sigma_Z,// this is Sigma (cov) not Omega (percision)
                          const arma::mat & y,
                          int k, int p, int n, 
                          int ns, int m, double emax // ars parameters
);


void update_Z_helper_Pois_reg(arma::mat & Z_curr, // persumably large, thus will not copy
                              const arma::mat & data, 
                              const arma::mat & design,
                              const arma::vec & mu_curr,
                              const arma::mat & beta_curr,
                              const arma::mat & Omega_curr,
                              int k, int p, int n, 
                              int ns, int m, double emax // ars parameters
);

void update_Z_helper_Pois_CAR(arma::mat & Z_curr, // persumably large, thus will not copy
                              const arma::mat & data, 
                              const arma::mat & design,
                              const arma::vec & mu_curr,
                              const arma::mat & beta_curr,
                              const arma::mat & Omega_curr,
                              int k, int p, int n, 
                              int ns, int m, double emax // ars parameters
);

void update_Z_helper_Pois_gra(arma::mat & Z_curr, // persumably large, thus will not copy
                              const arma::mat & data, 
                              const arma::vec & mu_curr,
                              const arma::mat & Omega_curr,
                              int k, int p, int n, 
                              int ns, int m, double emax // ars parameters
);

#endif
