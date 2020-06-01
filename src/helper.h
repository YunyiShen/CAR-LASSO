#ifndef HELPER_H
#define HELPER_H
arma::vec dLaplace_Cpp(const arma::vec & x, const double mu, const double lambda, bool take_log);
double rinvGau(double mu, const double & lambda);
arma::sp_mat getDesign_i_helper(const arma::rowvec & X_i,//row vector of that sample
                                int k);

#endif
