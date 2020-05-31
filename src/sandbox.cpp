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
double calculatee(int n){
  double res=0;
  Progress p(n, true);
  for(int i = 0 ; i < n ; ++i){
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return NA_REAL;
    } 
    res+=1/tgamma(i+1);
    p.increment();
  }
  return(res);
}