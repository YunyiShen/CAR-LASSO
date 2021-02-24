// transplanted from package bayesm

// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>

using namespace arma;
using namespace Rcpp;


List rwishart(double nu, arma::mat const& V){

// Wayne Taylor 4/7/2015

// Function to draw from Wishart (nu,V) and IW
 
// W ~ W(nu,V)
// E[W]=nuV

// WI=W^-1
// E[WI]=V^-1/(nu-m-1)
  
  // T has sqrt chisqs on diagonal and normals below diagonal
  int m = V.n_rows;
  arma::mat T = zeros(m,m);
  
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]); //rchisq returns a vectorized object, so using [0] allows for the conversion to double
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
  }}
  
  arma::mat C = trans(T)*chol(V);
  arma::mat CI = solve(trimatu(C),eye(m,m)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  // C is the upper triangular root of Wishart therefore, W=C'C
  // this is the LU decomposition Inv(W) = CICI' Note: this is
  // the UL decomp not LU!
  
  // W is Wishart draw, IW is W^-1
  
  return List::create(
    Named("W") = trans(C) * C,
    Named("IW") = CI * trans(CI),
    Named("C") = C,
    Named("CI") = CI);
}



// [[Rcpp::export]]
List rmultireg(arma::mat const& Y, arma::mat const& X, arma::mat const& Bbar, arma::mat const& A, double nu, arma::mat const& V) {

// Keunwoo Kim 09/09/2014

// Purpose: draw from posterior for Multivariate Regression Model with natural conjugate prior

// Arguments:
//  Y is n x m matrix
//  X is n x k
//  Bbar is the prior mean of regression coefficients  (k x m)
//  A is prior precision matrix
//  nu, V are parameters for prior on Sigma

// Output: list of B, Sigma draws of matrix of coefficients and Sigma matrix
 
// Model: 
//  Y=XB+U  cov(u_i) = Sigma
//  B is k x m matrix of coefficients

// Prior:  
//  beta|Sigma  ~ N(betabar,Sigma (x) A^-1)
//  betabar=vec(Bbar)
//  beta = vec(B) 
//  Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)

  int n = Y.n_rows;
  int m = Y.n_cols;
  int k = X.n_cols;
  
  //first draw Sigma
  arma::mat RA = chol(A);
  arma::mat W = join_cols(X, RA); //analogous to rbind() in R
  arma::mat Z = join_cols(Y, RA*Bbar);
  // note:  Y,X,A,Bbar must be matrices!
  arma::mat IR = solve(trimatu(chol(trans(W)*W)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  // W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
  arma::mat Btilde = (IR*trans(IR)) * (trans(W)*Z);
  // IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
  arma::mat E = Z-W*Btilde;
  arma::mat S = trans(E)*E;
  // E'E
  
  // compute the inverse of V+S
  arma::mat ucholinv = solve(trimatu(chol(V+S)), eye(m,m));
  arma::mat VSinv = ucholinv*trans(ucholinv);
  
  List rwout = rwishart(nu+n, VSinv);
  
  // now draw B given Sigma
  //   note beta ~ N(vec(Btilde),Sigma (x) Covxxa)
  //       Cov=(X'X + A)^-1  = IR t(IR)  
  //       Sigma=CICI'    
  //       therefore, cov(beta)= Omega = CICI' (x) IR IR' = (CI (x) IR) (CI (x) IR)'
  //  so to draw beta we do beta= vec(Btilde) +(CI (x) IR)vec(Z_mk)  
  //  		Z_mk is m x k matrix of N(0,1)
  //	since vec(ABC) = (C' (x) A)vec(B), we have 
  //		B = Btilde + IR Z_mk CI'

  arma::mat CI = rwout["CI"]; //there is no need to use as<arma::mat>(rwout["CI"]) since CI is being initiated as a arma::mat in the same line
  arma::mat draw = arma::mat(rnorm(k*m));
  draw.reshape(k,m);
  arma::mat B = Btilde + IR*draw*trans(CI);
    
  return List::create(
      Named("B") = B, 
      Named("Sigma") = rwout["IW"]);
}
