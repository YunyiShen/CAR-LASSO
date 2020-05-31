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
double intmin(int n){
  arma::vec res(n,fill::zeros);
  res += NA_REAL;
  
  if(R_IsNA(res(0))) return(1);
  return(0);
  
}