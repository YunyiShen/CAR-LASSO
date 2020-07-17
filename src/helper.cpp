// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;

// tested 20200601
arma::vec dLaplace_Cpp(const arma::vec & x, const double mu, const double lambda, bool take_log){
  arma::vec logd = - lambda * abs(x-mu) + log(lambda/2);
  return(take_log ? logd : exp(logd));
}


// inverse Gaussian random variable, tested already 20200528
/*
 * Mitchael,J.R., Schucany, W.R. and Haas, R.W. (1976). Generating
 random roots from variates using transformations with multiple roots.
 American Statistician. 30-2. 88-91.
 */	

// [[Rcpp::export]]
double rinvGau(double mu, double  lambda){
  mu = mu < 1e-12 ? 1e-12 : mu; // truncate mu
  lambda = lambda < 1e-12 ? 1e-12 : lambda;
  double b = 0.5 * mu / lambda;
  double a = mu * b;
  double c = 4.0 * mu * lambda;
  double d = mu * mu;
  double res;
  double u = R::runif(0,1);
  double chisqr = R::rchisq(1);
  double x = mu + a * chisqr - b * sqrt( c * chisqr + d * chisqr * chisqr); // solve the smaller root
  res = (u < (mu / (mu + x))) ? x : d/x; // accept the small one with prob mu/(mu+x), otherwise the larger one
  return(res);
}


// inverse Gaussian random variable, tested already 20200528
/*
 * Mitchael,J.R., Schucany, W.R. and Haas, R.W. (1976). Generating
 random roots from variates using transformations with multiple roots.
 American Statistician. 30-2. 88-91.
 */	
double rinvGau_full(const double & mu, const double & lambda){
  double b = 0.5 * mu / lambda;
  double a = mu * b;
  double c = 4.0 * mu * lambda;
  double d = mu * mu;
  double res;
  
  if (mu<=0 || lambda<=0) {
    return(NA_REAL);
  }
  
  double u = R::runif(0,1);
  double chisqr = R::rchisq(1);
  
  double x = mu + a * chisqr - b * sqrt( c * chisqr + d * chisqr * chisqr); // solve the smaller root
  res = (u < (mu / (mu + x))) ? x : d/x; // accept the small one with prob mu/(mu+x), otherwise the larger one
  return(res);
}


// get convinient X_i: design matrix for each sample, this will be blockwise diagnol, to be compactable for vector version of betas
// tested 20200528
arma::sp_mat getDesign_i_helper(const arma::rowvec & X_i,//row vector of that sample
                                int k){ // number of nodes
  int p = X_i.n_elem; // number of predictors
  arma::sp_mat Design(k,k * p);
  for(int j = 0 ; j < k ; ++j){
    Design(j,arma::span( j*p , (j+1) * p - 1 )) = X_i ; 
  }
  return(Design);
}

// resize a List
List resize( const List& x, int n ){
  int oldsize = x.size() ;
  List y(n) ;
  for( int i=0; i<oldsize; i++) y(i) = x(i) ;
  return y ;
}


// Inner function to simulate random uniforms in a matrix:
arma::mat RandMat(int nrow, int ncol)
{
  arma::mat Res = arma::randu(nrow,ncol);
  return(Res);
}


arma::sp_mat block_diag(const arma::mat & A, int d, int n){
  arma::sp_mat D(n*d,n*d);
  for(int i = 0 ; i < n ; ++i){
    D.submat(i*d, i*d, (i+1)*d-1, (i+1)*d-1) = A;
  }
  return(D);
}

// [[Rcpp::export]]
double stein_loss_cpp(const arma::mat & Omega, const arma::mat & Omega_hat){
  arma::mat prod_hat = solve(Omega_hat,Omega);
  int p = Omega.n_cols;
  double res = trace(prod_hat)-log(det(prod_hat))-p;
  return(res);
}
