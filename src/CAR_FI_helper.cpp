/*
 *This part was intend to calculate the FI matrix for CAR model, 
 *for experimental design aand active learning of such model
*/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// For coders: ALL THE INDECES START FROM 0!!

// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\beta_{ij}}\partial{beta}
//   result was a p by k matrix of above partial derivitives, useful in FI
//   need to be vectorized when constructing the FI
arma::mat sec_dev_betaij_beta(const arma::mat & design, // a row vector of design
                              const arma::mat & Sigma, // the cov mat
                              int i, int j){
    arma::mat XtX = design.t() * design;
    return(-XtX.col(i) * Sigma.row(j));
}

// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\beta_{ij}}\partial{mu}
//   result was a k by 1 matrix of above partial derivitives, useful in FI
arma::mat sec_dev_betaij_mu(const arma::mat & design, 
                              const arma::mat & Sigma,
                              int i, int j){
    return(-Sigma.col(j) * design(0,i));
}

// There is no need for partial l^2 / partial mu partial mu since it is -Sigma

// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\beta_{ij}}\partial{Omega(upper.tri)}
//   result was a k by 1 matrix of above partial derivitives, useful in FI
arma::mat sec_dev_betaij_Omega(const arma::mat & design, 
                               const arma::mat & Sigma,
                               const arma::mat & beta,
                               const arma::vec & mu,
                               int i, int j){

    arma::mat XtX = design.t() * design;
    arma::mat temp1,temp;
    temp1 = Sigma * mu * design; 
    temp = temp1.col(i) * Sigma.row(j) + Sigma * beta.t() * XtX.col(i) * Sigma.row(j);

    // get only the upper tri by chain rule
    temp += temp.t();
    temp.diag() /= 2;
    return(temp(trimatu_ind(size(temp),1)));
}


// This function calculate the second derivitave: 
//   \partial{l}^2/\partial{\mu_{j}}\partial{Omega(upper.tri)}
//   result was a k by 1 matrix of above partial derivitives, useful in FI
arma::mat sec_dev_betaij_Omega(const arma::mat & design, 
                               const arma::mat & Sigma,
                               const arma::mat & beta,
                               const arma::vec & mu,
                               int j){

    arma::mat temp;

    temp = Sigma.col(j) * (design * beta + mu.t()) * Sigma;

    // get only the upper tri by chain rule
    temp += temp.t();
    temp.diag() /= 2;
    return(temp(trimatu_ind(size(temp),1)));
}