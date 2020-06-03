// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
static const double pi = 3.141592653589793238462643383280;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "helper.h"

// Auto-normal

/* Starting here is SAR-LASSO and convert it to CAR if needed, see Ver Hoef et al. 2018
 * The Reason for using SAR rather than CAR is that SAR asks for a ralative easier to check condition 
 *   condition on B (I-B) non-singular, while CAR asks for C to have positive eigen value
 *   We can use independent Laplace prior in SAR to get similar result as LASSO.
 *   
 * This approach will allow us to do regressive LASSO and Graphical LASSO at the same time,
 *   while regressive LASSO ask known cov matrix and Graphical LASSO did not do anything to 
 *   expectation. 
 *   
 * A special Laplace prior was set here, which was related to the SAR noise sigma, while even we assume
 *   iid noise, SAR still can approimate any arbitrary covariance matrix, by Laplace prior to B matrix
 *   we will still get sparse solution to covariance matrix, which was the goal of Graphical LASSO
 *   
 * Method of fitting was Gibbs sampler, derived follow Park and Casella 2008 and Hoef et al. 2018
 * 
 */


// generate the B function from its vector version, we will not save diagonal since they are alwasy 0
// tested 20200528 // row first save

// [[Rcpp::export]]
arma::mat getBmat_helper(const arma::vec &B, int k){
  arma::mat Bmat(k,k,fill::zeros);
  arma::vec curr_col(k-1,fill::zeros);
  Bmat(arma::span(1,k-1),0) = B(arma::span(0,k-2)); //first column
  Bmat(arma::span(0,k-2),k-1) = B(arma::span((k-1) * (k-1),k * (k-1)-1)); // last column
  
  for(int i = 1 ; i < k-1 ; ++i){ // columns in between
    curr_col = B( arma::span(i * (k-1), (i+1)*(k-1)-1 ) );
    Bmat(arma::span(0,i-1),i) = curr_col(arma::span(0,i-1));
    Bmat(arma::span(i+1,k-1),i) = curr_col(arma::span(i,k-2));
  }
  return(Bmat.t());
}



// [[Rcpp::export]]
double Loglik_SAR_Cpp(const arma::vec & Z,
                      const arma::vec & mu,
                      const arma::vec & sigma,// sd here
                      const arma::mat B){
  int N = Z.n_elem;
  arma::mat I_minus_B = - B;
  I_minus_B.diag() += 1; // I - B
  //Rcout << I_minus_B << endl;
  arma::mat Q(N,N,fill::zeros);
  Q.diag() = 1 / sigma % sigma; // Var for noise
  
  // calculating normalizing constant
  // The cov matrix of SAR is : S = (I-B)^{-1} \Sigma (I^T-B^T)^{-1}
  // Thus log(det(S)) = log( det(\Sigma) ) - 2 log (det(I-B))
  double det_I_minus_B = det(I_minus_B);
  
  double logC = - (N/2) * log (2 * pi) - // 2pi
    sum(log(sigma)) + // Sigma^2s
    .5 * log((det_I_minus_B)*(det_I_minus_B)); // I-B part
    //Rcout << (det_I_minus_B)*(det_I_minus_B) << endl;
    arma::vec Y = Z-mu;
    
    
    double logH = as_scalar( - 0.5 * Y.t() * I_minus_B.t() * Q *  I_minus_B * Y); // Hamiltonian
    return(logC+logH);
}

// tested 20200602
// [[Rcpp::export]]
arma::mat SAR_to_Graph(const arma::vec & Bvec,const double & sigma){
  int k = floor( 0.5*( 1 + sqrt(1 + 4 * Bvec.n_elem)) );
  arma::mat B = getBmat_helper(Bvec, k);
  arma::mat Sigma = eye(size(B));
  Sigma *= 1/(sigma*sigma);
  arma::mat I_minus_B = -B;
  I_minus_B.diag() += 1;
  arma::mat full_percision = I_minus_B.t() * Sigma * I_minus_B; //cov matrix of this Sar model
  return(full_percision);
}


// [[Rcpp::export]]
Rcpp::List Sigma_to_CAR_Cpp(const arma::mat & Sigma){
  
  arma::mat Q = inv(Sigma); // percision matrix
  arma::mat R = -Sigma; // R matrix, see Ver Hoef et al. 2018
  R.diag() *= 0;
  arma::vec M = 1/Q.diag();
  arma::mat invD(size(Q));
  invD.diag() = M ;
  return(Rcpp::List::create(
      Rcpp::Named("C") = invD * R, // C matrix per Ver Hoef et al. 2018
      Rcpp::Named("M") = M // diag of M matrix per Ver Hoef et al. 2018
  ));
}


// tested 20200529
// [[Rcpp::export]]
arma::mat Sample_SAR_Cpp(const int Nsample,
                         const arma::vec & mu,
                         const arma::vec & sigma, // make sure this is standard dev
                         const arma::mat B){
  int N = mu.n_elem;
  arma::mat I_minus_B = - B;
  I_minus_B.diag() += 1; // I - B
  arma::mat Sigma(N,N,fill::zeros);
  Sigma.diag() = sigma%sigma; // Var for noise
  arma::mat try_solve;
  inv(try_solve , I_minus_B);// This will be all false if I_minus_B is singular
  if(!try_solve(0,0)){
    return(try_solve); // if I-B not inversible, return 0s
  }
  
  arma::mat full_Sigma = try_solve * Sigma * try_solve.t(); //cov matrix of this Sar model
  arma::mat Res;
  Res = mvnrnd(mu, full_Sigma, Nsample ); // this will give column vectors
  return(Res);
}


// [[Rcpp::export]]
arma::mat Sample_SAR_from_Design_Cpp(const arma::mat & Design,
                                     const arma::mat & beta,
                                     const arma::vec & mu,
                                     const arma::vec & sigma,
                                     const arma::mat & B){
  int n = Design.n_rows;
  int k = B.n_rows;
  arma::mat Res(k,n,fill::zeros);
  
  arma::mat thr = Design * beta;
  thr = thr.t();
  thr.each_col() += mu;
  
  for(int i = 0 ; i < n ; ++i){
    Res.col(i) = Sample_SAR_Cpp(1,thr.col(i),sigma,B);
  }
  return(Res);
}

