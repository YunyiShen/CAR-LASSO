// [[Rcpp::depends(RcppArmadillo)]]
#include <tgmath.h>
#include <RcppArmadillo.h> 
#include <iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;



// [[Rcpp::export]]
double gig_mode(double lambda, double omega)
{
  if (lambda >= 1.)
    // mode of fgig(x) //
    return (sqrt((lambda-1.)*(lambda-1.) + omega*omega)+(lambda-1.))/omega;
  else
    // 0 <= lambda < 1: use mode of f(1/x) //
    return omega / (sqrt((1.-lambda)*(1.-lambda) + omega*omega)+(1.-lambda));
} 


// [[Rcpp::export]]
double rgig_ROU_noshift(double lambda, 
                        double lambda_old, 
                        double omega, 
                        double alpha){
  double xm, nc;     // location of mode; c=log(f(xm)) normalization constant //
  double ym, um;     // location of maximum of x*sqrt(f(x)); umax of MBR //
  double s, t;       // auxiliary variables //
  double U, V, X;    // random variables //
  double res=0;
  //int i;             // loop variable (number of generated random variables) //
  int count = 0;     // counter for total number of iterations //
  
  t = 0.5 * (lambda-1.0);
  s = 0.25 * omega;
  
  // mode = location of maximum of sqrt(f(x)) //
  xm = gig_mode(lambda, omega);
  
  // normalization constant: c = log(sqrt(f(xm))) //
  nc = t*log(xm) - s*(xm + 1.0/xm);
  
  // location of maximum of x*sqrt(f(x)):           //
  // we need the positive root of                   //
  //    omega/2*y^2 - (lambda+1)*y - omega/2 = 0    //
  ym = ((lambda+1.0) + sqrt((lambda+1.0)*(lambda+1.0) + omega*omega))/omega;
  
  // boundaries of minmal bounding rectangle:                   //
  // we us the "normalized" density f(x) / f(xm). hence         //
  // upper boundary: vmax = 1.                                  //
  // left hand boundary: umin = 0.                              //
  // right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) //
  um = exp(0.5*(lambda+1.0)*log(ym) - s*(ym + 1.0/ym) - nc);
  

  do {
    ++count;
    U = um * R::runif(0,1);        // U(0,umax) //
    V = R::runif(0,1);             // U(0,vmax) //
    X = U/V;
  }                              // Acceptance/Rejection //
  while (((log(V)) > (t*log(X) - s*(X + 1.0/X) - nc)));
    
    // store random point //
    res = (lambda_old < 0.0) ? (alpha / X) : (alpha * X);
  
  return(res);
} 








