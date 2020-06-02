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
// [[Rcpp::export]]
arma::sp_mat block_diag(const arma::mat A, int d, int n){
  arma::sp_mat D(n*d,n*d);
  for(int i = 0 ; i < n ; ++i){
    D.submat(i*d, i*d, (i+1)*d-1, (i+1)*d-1) = A;
  }
  return(D);
}

