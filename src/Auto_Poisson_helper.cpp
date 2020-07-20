// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"

/* Convention:
 *    single sample will be saved as col vectors
 */

// Auto-Poisson

arma::vec Auto_Poisson_Gibbs_Single_Cpp(const arma::sp_mat &graph, // graph
                                        const arma::vec &thresholds, // mean 
                                        const int & Winsorized, // right censor
                                        const int & nIter){  // number of Gibbs iters
  int N = graph.n_rows;
  int temp;
  arma::vec Res(N);
  
  for (int i = 0 ; i < N ; ++i){
    Res(i) = R::rpois(thresholds(i)); // initializing 
  }
  
  for(int i = 0 ; i < nIter ; ++i){ // outer loop, for number of samples
    for(int j = 0 ; j < N ; ++j){
      temp = R::rpois(exp(thresholds(j) + 0.5 * as_scalar(Res.t() * graph.col(j)))); // Gibbs scan
      Res(i) = temp > Winsorized ? Winsorized : temp; // right censored 
    }
  }
  return(Res);
} 


// [[Rcpp::export]]
arma::mat Auto_Poisson_Gibbs_Batch_Cpp(const int & Nsample, // number of samples
                                       const arma::sp_mat &graph, // graph
                                       const arma::vec &thresholds, // mean 
                                       const int & Winsorized, // right censor
                                       const int & nIter){  
  int N = graph.n_rows;
  arma::mat Res(N,Nsample);
  for(int i = 0 ; i < Nsample ; ++i ){
    Res.col(i) = Auto_Poisson_Gibbs_Single_Cpp(graph,thresholds,Winsorized,nIter);
  }
  return(Res);
}

// [[Rcpp::export]]
double Auto_Poisson_Pseudo_likelihood_Cpp(const arma::mat Sample,
                                          const arma::sp_mat & graph,
                                          const arma::vec & thresholds){
  int Nsample = Sample.n_cols; //samples
  int N = Sample.n_rows; // nodes
  double Res = 0;
  arma::mat Interactions = graph * Sample ; // interaction terms
  for(int i = 0 ; i < Nsample ; ++i){
    arma::vec Lambda = exp(thresholds + Interactions.col(i)); // lambda of Poisson of that repeat
    for(int j = 0 ; j < N ; ++j){
      Res += R::dpois(Sample(j,i),Lambda(j),true); // likelihood, use Poisson approximation 
    } 
  }
  
  return(Res);
}

// Q-function for auto poisson 
double Auto_Poisson_Q_Cpp(const arma::mat Sample,
                          const arma::sp_mat & graph,
                          const arma::vec & thresholds){
  // not sure if this line work sonce sum(...) will give a row vec
  arma::rowvec thr_contribution = thresholds.t() * Sample - sum( lgamma(Sample + 1),0); // non-interaction part of Besag's Q-function, see Augustin et al. 2006
  double Res = 0;
  int Nsample = Sample.n_cols;
  
  for(int i = 0 ; i < Nsample ; ++i){
    Res += thr_contribution(i) + .5 * as_scalar(Sample.col(i).t() * graph * Sample.col(i)); // Besag's Q function for (spatial) auto Poisson model 
  }
  return(Res);
}


// logLik of a normal response using composition as predictor
double logLik_Normal_Response(const arma::vec & data,
                              const arma::mat & Design,
                              const arma::vec & beta_fix,
                              const arma::mat & composition,
                              const arma::vec & beta_composition,
                              const double & sigma){
  int N = Design.n_rows;
  arma::vec mu = Design * beta_fix + composition.t() * beta_composition;
  double Res = 0;
  for(int i = 0 ; i < N ; ++i){
    Res += R::dnorm(data(i) , mu(i) , sigma,true);
  }
  return(Res);
}





/* TODO:
 * 1. Coupling from the past for one perfect sample from the auto poisson model
 * 2. Auto-NegBin model sampling, 
 * 3. derive the Q-function 
 */







