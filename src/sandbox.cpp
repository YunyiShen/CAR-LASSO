// [[Rcpp::depends(RcppArmadillo)]]
#include <tgmath.h>
#include <RcppArmadillo.h> 
#include <iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;


// generate the B function from its vector version, we will not save diagonal since they are alwasy 0
// tested 20200528 
arma::mat gettingBmat_helper(const arma::vec &B, int k){
  arma::mat Bmat(k,k,fill::zeros);
  arma::vec curr_col(k-1,fill::zeros);
  Bmat(span(1,k-1),0) = B(span(0,k-2)); //first column
  Bmat(span(0,k-2),k-1) = B(span((k-1) * (k-1),k * (k-1)-1)); // last column
  
  for(int i = 1 ; i < k-1 ; ++i){ // columns in between
    curr_col = B( span(i * (k-1), (i+1)*(k-1)-1 ) );
    Bmat(span(0,i-1),i) = curr_col(span(0,i-1));
    Bmat(span(i+1,k-1),i) = curr_col(span(i,k-2));
  }
  return(Bmat);
}

// get convinient X_i: design matrix for each sample
// tested 20200528
arma::sp_mat getDesign_i_helper(const arma::rowvec & X_i,//row vector of that sample
                                int k){ // number of nodes
  int p = X_i.n_elem; // number of predictors
  arma::sp_mat Design(k,k * p);
  for(int j = 0 ; j < k ; ++j){
    Design(j,span( j*p , (j+1) * p - 1 )) = X_i ; 
  }
  return(Design);
}




// [[Rcpp::export]]
arma::vec update_B_helper(const arma::mat & data,
                          const arma::mat & design,
                          const arma::vec & mu,
                          const arma::vec & beta,
                          const double & sigma2,
                          const arma::vec & eta2,
                          int k, int n){
  arma::mat Omega(k,k,fill::zeros);
  arma::vec Y_i_tilde;
  arma::vec B(k*(k-1));
  for(int i = 0 ; i < n ; ++i){
    Y_i_tilde = getDesign_i_helper(design.row(i),k) * beta - data.col(i) + mu;
    Omega += Y_i_tilde * Y_i_tilde.t(); // The 
  }
  
  arma::mat Omega_temp;
  arma::mat Sigma_temp;
  arma::mat mu_i;
  arma::uvec ind = linspace<uvec>(0, k-1 , k); // this is a index vector, for deleting the digonal related row and colunms in var cov mat
  arma::uvec ind_remain;
  for(int i = 0 ; i < k ; ++i){
    ind_remain = find(ind != i);
    Omega_temp = Omega.submat(ind_remain,ind_remain); // delete row and column related with the diagnol entry of B
    mu_i = Omega.col(i); // used to solve for expectation
    mu_i = mu_i(ind_remain);
    Omega_temp.diag() += 1/(eta2(span( i*(k-1) , (i+1)*(k-1) - 1 ))); // add LASSO shrik
    
    Sigma_temp = inv(Omega_temp); // calculate cov mat up to sigma^2 scale
    
    mu_i = Sigma_temp * mu_i; // sove for expectation
    
    Sigma_temp = sigma2 * Sigma_temp; // real cov mat
    
    
    B( span(i*(k-1) , (i+1)*(k-1) - 1) ) = mvnrnd(mu_i, Sigma_temp); // save in vector
  }
  
  return(B);
}






