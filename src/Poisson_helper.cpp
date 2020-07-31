// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
#include <progress.hpp>
#include <progress_bar.hpp>

/* This calculate the log Posterior (up to normalization) of latent normal variable Zij condition 
 * on Omega mu and all other Zs. 
 * 
 * Also calculate the derivitave, res(0) is the value while res(1) is the derivative
 * 
 * This function will be used in adaptive rejection sampling, since it is log concave
 */

arma::vec logPostZij_helper_Cpp(int i, int j,// which node
                                const arma::mat & Z_curr,
                                const arma::mat & mu_Z,
                                const arma::mat & Sigma_Z,// this is Sigma (cov) not Omega (percision)
                                const arma::mat & y,
                                int k, int p, int n){
  arma::vec res(2,fill::zeros);
  res(0) += y(i,j) * Z_curr(i,j) - exp(Z_curr(i,j)); // log posterior due to Poisson
  res(1) += y(i,j) - exp(Z_curr(i,j)); // d/dz logPost due to Poisson
  
  arma::uvec ind = linspace<uvec>(0,k-1,k);
  arma::uvec ind_noi = find(ind!=i);
  arma::uvec ind_i = find(ind==i);
  arma::uvec ind_j = find(ind==j);
  
  double Sigmabb = Sigma_Z(i,i);
  arma::mat Sigmac = Sigma_Z(ind_noi,ind_i);
  arma::mat Sigmaa = Sigma_Z(ind_noi,ind_noi);
  
  double mu_Zij = mu_Z(i) + as_scalar(Sigmac.t() * solve( Sigmaa ,trans(Z_curr(ind_j,ind_noi)-mu_Z(ind_j,ind_noi))));
  double sigma2_Zij = Sigmabb - as_scalar( Sigmac.t() * solve(Sigmaa,Sigmac));
  
  
  res(0) -= 0.5 * (Z_curr(i,j)-mu_Zij)*(Z_curr(i,j)-mu_Zij)/sigma2_Zij;
  res(1) -= (Z_curr(i,j)-mu_Zij)/sigma2_Zij;
  return(res);
}




