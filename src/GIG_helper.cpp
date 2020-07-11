// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

/*
 * Adopted form R package GIGrvg
 * 
 * 
 */


//***************************************************************************//
// Privat Functions                                                          //
//***************************************************************************//

double _gig_mode(double lambda, double omega)
  
  // Compute mode of GIG distribution.                                         //
  //                                                                           //
  // Parameters:                                                               //
  //   lambda .. parameter for distribution                                    //
  //   omega ... parameter for distribution                                    //
  //                                                                           //
  // Return:                                                                   //
  //   mode                                                                    //
  
{
  if (lambda >= 1.)
    // mode of fgig(x) //
    return (sqrt((lambda-1.)*(lambda-1.) + omega*omega)+(lambda-1.))/omega;
  else
    // 0 <= lambda < 1: use mode of f(1/x) //
    return omega / (sqrt((1.-lambda)*(1.-lambda) + omega*omega)+(1.-lambda));
} // end of _gig_mode() //



void _rgig_ROU_noshift (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
  
  // Tpye 1:                                                                   //
  // Ratio-of-uniforms without shift.                                          //
  //   Dagpunar (1988), Sect.~4.6.2                                            //
  //   Lehner (1989)                                                           //
  
{
  double xm, nc;     // location of mode; c=log(f(xm)) normalization constant //
  double ym, um;     // location of maximum of x*sqrt(f(x)); umax of MBR //
  double s, t;       // auxiliary variables //
  double U, V, X;    // random variables //
  
  int i;             // loop variable (number of generated random variables) //
  int count = 0;     // counter for total number of iterations //
  
  // -- Setup -------------------------------------------------------------- //
  
  // shortcuts //
  t = 0.5 * (lambda-1.);
  s = 0.25 * omega;
  
  // mode = location of maximum of sqrt(f(x)) //
  xm = _gig_mode(lambda, omega);
  
  // normalization constant: c = log(sqrt(f(xm))) //
  nc = t*log(xm) - s*(xm + 1./xm);
  
  // location of maximum of x*sqrt(f(x)):           //
  // we need the positive root of                   //
  //    omega/2*y^2 - (lambda+1)*y - omega/2 = 0    //
  ym = ((lambda+1.) + sqrt((lambda+1.)*(lambda+1.) + omega*omega))/omega;
  
  // boundaries of minmal bounding rectangle:                   //
  // we us the "normalized" density f(x) / f(xm). hence         //
  // upper boundary: vmax = 1.                                  //
  // left hand boundary: umin = 0.                              //
  // right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) //
  um = exp(0.5*(lambda+1.)*log(ym) - s*(ym + 1./ym) - nc);
  
  // -- Generate sample ---------------------------------------------------- //
  
  for (i=0; i<n; i++) {
    do {
      ++count;
      U = um * R::runif(0,1);        // U(0,umax) //
      V = R::runif(0,1);             // U(0,vmax) //
      X = U/V;
    }                              // Acceptance/Rejection //
    while (((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));
    
    // store random point //
    res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
  }
  
  // -- End ---------------------------------------------------------------- //
  
  return;
} // end of _rgig_ROU_noshift() //




void _rgig_newapproach1 (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
  
  // Type 4:                                                                   //
  // New approach, constant hat in log-concave part.                           //
  // Draw sample from GIG distribution.                                        //
  //                                                                           //
  // Case: 0 < lambda < 1, 0 < omega < 1                                       //
  //                                                                           //
  // Parameters:                                                               //
  //   n ....... sample size (positive integer)                                //
  //   lambda .. parameter for distribution                                    //
  //   omega ... parameter for distribution                                    //
  //                                                                           //
  // Return:                                                                   //
  //   random sample of size 'n'                                               //
  
{
  // parameters for hat function //
  double A[3], Atot;  // area below hat //
  double k0;          // maximum of PDF //
  double k1, k2;      // multiplicative constant //
  
  double xm;          // location of mode //
  double x0;          // splitting point T-concave / T-convex //
  double a;           // auxiliary variable //
  
  double U, V, X;     // random numbers //
  double hx;          // hat at X //
  
  int i;              // loop variable (number of generated random variables) //
  int count = 0;      // counter for total number of iterations //
  
  // -- Check arguments ---------------------------------------------------- //
  
  if (lambda >= 1. || omega >1.)
    stop("invalid parameters");
  
  // -- Setup -------------------------------------------------------------- //
  
  // mode = location of maximum of sqrt(f(x)) //
  xm = _gig_mode(lambda, omega);
  
  // splitting point //
  x0 = omega/(1.-lambda);
  
  // domain [0, x_0] //
  k0 = exp((lambda-1.)*log(xm) - 0.5*omega*(xm + 1./xm));     // = f(xm) //
  A[0] = k0 * x0;
  
  // domain [x_0, Infinity] //
  if (x0 >= 2./omega) {
    k1 = 0.;
    A[1] = 0.;
    k2 = pow(x0, lambda-1.);
    A[2] = k2 * 2. * exp(-omega*x0/2.)/omega;
  }
  
  else {
    // domain [x_0, 2/omega] //
    k1 = exp(-omega);
    A[1] = (lambda == 0.) 
      ? k1 * log(2./(omega*omega))
        : k1 / lambda * ( pow(2./omega, lambda) - pow(x0, lambda) );
    
    // domain [2/omega, Infinity] //
    k2 = pow(2/omega, lambda-1.);
    A[2] = k2 * 2 * exp(-1.)/omega;
  }
  
  // total area //
  Atot = A[0] + A[1] + A[2];
  
  // -- Generate sample ---------------------------------------------------- //
  
  for (i=0; i<n; i++) {
    do {
      ++count;
      
      // get uniform random number //
      V = Atot * R::runif(0,1);
      
      do {
        
        // domain [0, x_0] //
        if (V <= A[0]) {
          X = x0 * V / A[0];
          hx = k0;
          break;
        }
        
        // domain [x_0, 2/omega] //
        V -= A[0];
        if (V <= A[1]) {
          if (lambda == 0.) {
            X = omega * exp(exp(omega)*V);
            hx = k1 / X;
          }
          else {
            X = pow(pow(x0, lambda) + (lambda / k1 * V), 1./lambda);
            hx = k1 * pow(X, lambda-1.);
          }
          break;
        }
        
        // domain [max(x0,2/omega), Infinity] //
        V -= A[1];
        a = (x0 > 2./omega) ? x0 : 2./omega;
        X = -2./omega * log(exp(-omega/2. * a) - omega/(2.*k2) * V);
        hx = k2 * exp(-omega/2. * X);
        break;
        
      } while(0);
      
      // accept or reject //
      U = unif_rand() * hx;
      
      if (log(U) <= (lambda-1.) * log(X) - omega/2. * (X+1./X)) {
        // store random point //
        res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
        break;
      }
    } while(1);
    
  }
  
  // -- End ---------------------------------------------------------------- //
  
  return;
} 
// end of _rgig_newapproach1() //



void _rgig_ROU_shift_alt (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
  
  // Type 8:                                                                   //
  // Ratio-of-uniforms with shift by 'mode', alternative implementation.       //
  //   Dagpunar (1989)                                                         //
  //   Lehner (1989)                                                           //
  
  
{
  double xm, nc;     // location of mode; c=log(f(xm)) normalization constant //
  double s, t;       // auxiliary variables //
  double U, V, X;    // random variables //
  
  int i;             // loop variable (number of generated random variables) //
  int count = 0;     // counter for total number of iterations //
  
  double a, b, c;    // coefficent of cubic //
  double p, q;       // coefficents of depressed cubic //
  double fi, fak;    // auxiliary results for Cardano's rule //
  
  
  double y1, y2;     // roots of (1/x)*sqrt(f((1/x)+m)) //
  
  double uplus, uminus;  // maximum and minimum of x*sqrt(f(x+m)) //
  
  t = 0.5 * (lambda-1.);
  s = 0.25 * omega;
  
  // mode = location of maximum of sqrt(f(x)) //
  xm = _gig_mode(lambda, omega);
  
  // normalization constant: c = log(sqrt(f(xm))) //
  nc = t*log(xm) - s*(xm + 1./xm);
  
  // location of minimum and maximum of (1/x)*sqrt(f(1/x+m)):  //
  
  // compute coeffients of cubic equation y^3+a*y^2+b*y+c=0 //
  a = -(2.*(lambda+1.)/omega + xm);       // < 0 //
  b = (2.*(lambda-1.)*xm/omega - 1.);
  c = xm;
  
  // we need the roots in (0,xm) and (xm,inf) //
  
  // substitute y=z-a/3 for depressed cubic equation z^3+p*z+q=0 //
  p = b - a*a/3.;
  q = (2.*a*a*a)/27. - (a*b)/3. + c;
  
  // use Cardano's rule //
  fi = acos(-q/(2.*sqrt(-(p*p*p)/27.)));
  fak = 2.*sqrt(-p/3.);
  y1 = fak * cos(fi/3.) - a/3.;
  y2 = fak * cos(fi/3. + 4./3.*M_PI) - a/3.;
  
  // boundaries of minmal bounding rectangle:                  //
  // we us the "normalized" density f(x) / f(xm). hence        //
  // upper boundary: vmax = 1.                                 //
  // left hand boundary: uminus = (y2-xm) * sqrt(f(y2)) / sqrt(f(xm)) //
  // right hand boundary: uplus = (y1-xm) * sqrt(f(y1)) / sqrt(f(xm)) //
  uplus  = (y1-xm) * exp(t*log(y1) - s*(y1 + 1./y1) - nc);
  uminus = (y2-xm) * exp(t*log(y2) - s*(y2 + 1./y2) - nc);
  
  // -- Generate sample ---------------------------------------------------- //
  
  for (i=0; i<n; i++) {
    do {
      ++count;
      U = uminus + R::runif(0,1) * (uplus - uminus);    // U(u-,u+)  //
      V = R::runif(0,1);                                // U(0,vmax) //
      X = U/V + xm;
    }                                         // Acceptance/Rejection //
    while ((X <= 0.) || ((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));
    
    // store random point //
    res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
  }
  return;
} 


// Draw sample from GIG distribution.                                        //
// Wrapper for do_rgig() with GetRNGstate() ... PutRNGstate()                //
//                                                                           //
// Parameters:                                                               //
//   n ....... sample size (positive integer)                                //
//   lambda .. parameter for distribution                                    //
//   chi   ... parameter for distribution                                    //
//   psi   ... parameter for distribution                                    //
//                                                                           //
// Return:                                                                   //
//   random sample of size 'n'                                               //



// [[Rcpp::export]]
double rgig(double lambda, double chi, double psi){
  double omega, alpha;     // parameters of standard distribution //
  double res;
  
  // allocate array for random sample //
  
  if (chi < DOUBLE_EPS*10.0) { 
    // special cases which are basically Gamma and Inverse Gamma distribution //
    if (lambda > 0.0) {
      res = R::rgamma(lambda, 2.0/psi); 
    }
    else {
      res = 1.0/R::rgamma(-lambda, 2.0/psi); 
    }    
  }
  
  else if (psi < DOUBLE_EPS*10.0) {
    // special cases which are basically Gamma and Inverse Gamma distribution //
    if (lambda > 0.0) {
      res = 1.0/R::rgamma(lambda, 2.0/chi); 
    }
    else {
      res = R::rgamma(-lambda, 2.0/chi); 
    }    
    
  }
  
  else {
    double lambda_old = lambda;
    if (lambda < 0.) lambda = -lambda;
    alpha = sqrt(chi/psi);
    omega = sqrt(psi*chi);
    
    // run generator //
    do {
      if (lambda > 2. || omega > 3.) {
        // Ratio-of-uniforms with shift by 'mode', alternative implementation //
        _rgig_ROU_shift_alt(&res, 1, lambda, lambda_old, omega, alpha);
        break;
      }
      
      if (lambda >= 1.-2.25*omega*omega || omega > 0.2) {
        // Ratio-of-uniforms without shift //
        _rgig_ROU_noshift(&res, 1, lambda, lambda_old, omega, alpha);
        break;
      }
      
      if (lambda >= 0. && omega > 0.) {
        // New approach, constant hat in log-concave part. //
        _rgig_newapproach1(&res, 1, lambda, lambda_old, omega, alpha);
        break;
      }
      
      // else //
      stop("parameters must satisfy lambda>=0 and omega>0.");
      
    } while (0);
  }
  
  // return result //
  return res;
  
} 
  
