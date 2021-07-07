#ifndef CAR_ALASSO_RANDEFF_HELPER_H
#define CAR_ALASSO_RANDEFF_HELPER_H

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
#include "GIG_helper.h"
#include "CAR_LASSO_helper.h"
#include "CAR_LASSO_randeff_helper.h"

void update_car_randeff_Omega_adp_helper(arma::mat & Omega,
                             const arma::mat & data,
                             const arma::mat & design,
                             const arma::mat & design_r,
                             const arma::mat & nu, 
                             const arma::vec & mu,
                             const arma::mat & beta,
                             const arma::vec & lambda_curr,// different lambda for different entries
                             const arma::vec & lambda_diag, 
                             int k, int p,int n);



#endif