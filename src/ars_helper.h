#ifndef ARS_HELPER_H
#define ARS_HELPER_H
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;
//#include "Error.h"
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <iostream>
#include <climits>
#include <cmath>


void initial_(const int *ns, const int *m, const double *emax, const double *x, const double *hx, const double *hpx, 
                const int *lb, double *xlb,  const int *ub,      double *xub,     int *ifault,      int *iwv,           double *rwv);

void splhull_(double *u2,          int *ipt,       int *ilow,    const int *lb,  double *xlb, 
              double *hulb,        double *huzmax, double *alcu, double *x,      double *hx, 
              double *hpx,         double *z__,    double *huz,  double *scum,   const double *eps, 
              const double *emax,  double *beta,   int *i__,     int *j);


void splhull_(double *u2,          int *ipt,       int *ilow,    const int *lb,  double *xlb, 
              double *hulb,        double *huzmax, double *alcu, double *x,      double *hx, 
              double *hpx,         double *z__,    double *huz,  double *scum,   const double *eps, 
              const double *emax,  double *beta,   int *i__,     int *j);

void intersection_(double *x1,  double *y1, double *yp1, double *x2,        double *y2, 
                   double *yp2, double *z1, double *hz1, const double *eps, int *ifault);


void update_(int *n,         int *ilow,          int *ihigh,          int *ipt,       double *scum, 
             double *cu,     double *x,          double *hx,          double *hpx,    double *z__, 
             double *huz,    double *huzmax,     const double *emax, 
             const int *lb,  double *xlb,        double *hulb,        const int *ub,  double *xub,   double *huub, 
             int *ifault,    const double *eps,  double *alcu);

inline double 
expon_(const double *x, const double *emax)            /* exponential without underflow */
{
  return (*x < -(*emax)) ? 0.0 : exp(*x);
} 

#endif