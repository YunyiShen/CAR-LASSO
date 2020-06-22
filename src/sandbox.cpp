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
#include "helper.h"
#include "SRG_LASSO_helper.h"
#include "Truncated_Normal_helper.h"


void update(double * sup){
  --sup;
  *(sup+1) = 10000;
  return;
}
// [[Rcpp::export]]
arma::vec testtt(arma::vec & foo, int i){
  update(&foo(i));
  return(foo);
}







