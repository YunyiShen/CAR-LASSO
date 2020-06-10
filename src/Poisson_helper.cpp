// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h> // to use sparse matrix
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
  
  arma::uvec ind = linspace(0,k-1,k);
  arma::uvec ind_noi = find(ind!=i);
  
  Sigmabb = Sigma_Z(i,i);
  Sigmac = Sigma_Z(ind_noi,i);
  Sigmaa = Sigma_Z(ind_noi,ind_noi);
  
  mu_Zij = mu_Z(i) + t(Sigmac) * solve( Sigmaa ,t(Z_curr(j,ind_noi)-mu_Z(j,ind_noi)));
  sigma2_Zij = Sigmabb - t(Sigmac) * solve(Sigmaa,Sigmac);
  
  mu_Zij = as_scalar(mu_Zij);
  sigma2_Zij = as_scalar(sigma_Zij);
  
  res(0) -= 0.5 * (Z_curr(i,j)-mu_Zij)*(Z_curr(i,j)-mu_Zij)/sigma2_Zij;
  res(1) -= (Z_curr(i,j)-mu_Zij)/sigma2_Zij;
  
}