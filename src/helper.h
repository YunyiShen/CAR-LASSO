#ifndef HELPER_H
#define HELPER_H
// Laplace distribution density 
static const double pi = 3.141592653589793238462643383280;

arma::vec dLaplace_Cpp(const arma::vec & x, 
                       const double mu, 
                       const double lambda, 
                       bool take_log);

// sample from inverse Gaussian 
double rinvGau(double mu, double  lambda);

// get blockwise diagnoal design matrix
arma::sp_mat getDesign_i_helper(const arma::rowvec & X_i,//row vector of that sample
                                int k);

// get a random matrix
arma::mat RandMat(int nrow, int ncol);

// resize a List
Rcpp::List resize( const Rcpp::List& x, int n );

// blockwise diagnonal matrix from block
arma::sp_mat block_diag(const arma::mat & A, int d,int n);

#endif
