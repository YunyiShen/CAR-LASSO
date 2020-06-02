#ifndef HELPER_H
#define HELPER_H
arma::vec dLaplace_Cpp(const arma::vec & x, const double mu, const double lambda, bool take_log);
double rinvGau(double mu, double lambda);
arma::sp_mat getDesign_i_helper(const arma::rowvec & X_i,//row vector of that sample
                                int k);
arma::mat RandMat(int nrow, int ncol);
Rcpp::List resize( const Rcpp::List& x, int n );
arma::sp_mat block_diag(const arma::mat A, int d,int n);

#endif
