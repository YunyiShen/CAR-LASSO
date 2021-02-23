#ifndef ARS_MULTINOMIAL_HELPER_H
#define ARS_MULTINOMIAL_HELPER_H
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

void sample_multi_(int *iwv,      double *rwv, 
             int *ifault, int l, int w,// which node
             arma::mat & Z_curr,
             const arma::mat & mu_Z,
             const arma::mat & Sigma_Z,// this is Sigma (cov) not Omega (percision)
             const arma::mat & y,
             int k, int p, int n_sample);


void update_Z_helper_multinomial(arma::mat & Z_curr,
                          const arma::mat & mu_Z,
                          const arma::mat & Sigma_Z,// this is Sigma (cov) not Omega (percision)
                          const arma::mat & y,
                          int k, int p, int n, 
                          int ns, int m, double emax // ars parameters
);


void update_Z_helper_multinomial_SRG(arma::mat & Z_curr, // persumably large, thus will not copy
                              const arma::mat & data, 
                              const arma::mat & design,
                              const arma::vec & mu_curr,
                              const arma::mat & beta_curr,
                              const arma::mat & Omega_curr,
                              int k, int p, int n, 
                              int ns, int m, double emax // ars parameters
);

void update_Z_helper_multinomial_CAR(arma::mat & Z_curr, // persumably large, thus will not copy
                              const arma::mat & data, 
                              const arma::mat & design,
                              const arma::vec & mu_curr,
                              const arma::mat & beta_curr,
                              const arma::mat & Omega_curr,
                              int k, int p, int n, 
                              int ns, int m, double emax // ars parameters
);

void update_Z_helper_multinomial_gra(arma::mat & Z_curr, // persumably large, thus will not copy
                              const arma::mat & data, 
                              const arma::vec & mu_curr,
                              const arma::mat & Omega_curr,
                              int k, int p, int n, 
                              int ns, int m, double emax // ars parameters
);

// void update_Z_helper_multinomial_para(arma::mat &Z_curr,
//                                  const arma::mat &mu_Z,
//                                  const arma::mat &Sigma_Z, // this is Sigma (cov) not Omega (percision)
//                                  const arma::mat &y,
//                                  int k, int p, int n,
//                                  int ns, int m, double emax // ars parameters
// );

void update_Z_helper_multinomial_CAR_randeff(arma::mat &Z_curr, // persumably large, thus will not copy
                                     const arma::mat &data,
                                     const arma::mat &design,
                                     const arma::mat &design_r,
                                     const arma::vec &mu_curr,
                                     const arma::mat &beta_curr,
                                     const arma::mat &nu_curr,
                                     const arma::mat &Omega_curr,
                                     int k, int p, int n,
                                     int ns, int m, double emax // ars parameters
);

#endif
