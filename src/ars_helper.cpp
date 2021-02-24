/*** ars.cpp ***/
//
//    AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//    Modificado por: Paulino Perez Rodriguez
//    Fecha: 17/02/07
//    PURPOSE: Adaptive rejection sampling
//
/* ********************************************************************************* */
// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>
// #include <tgmath.h>
// using namespace Rcpp;
// using namespace arma;
// //#include "Error.h"
// #include <R.h>
// #include <Rmath.h>
// #include <Rdefines.h>
// #include <iostream>
// #include <climits>
// #include <cmath>
#include "ars_helper.h"

// Subroutine initial_:
// ====================
/* This subroutine takes as input the number of starting values m  */
/*  and the starting values x(i), hx(i), hpx(i)  i=1,m             */
/* As output we have pointer ipt along with ilow and ihigh and the lower */
/* and upper hulls defined  by z, hz, scum, cu, hulb, huub stored in working */
/* vectors iwv and rwv */
/* Ifault detects wrong starting points or non-concavity */

/* DESCRIPTION OF PARAMETERS and place of storage                                */
/*                                                                               */
/*     ns   iwv(3) : maximum number of points defining the hulls                 */
/*     m    iwv(4) : number of starting points                                   */
/*     emax rwv(3) : large value for which it is possible to compute             */
/*                   an exponential, eps=exp(-emax) is taken as a small          */
/*                   value used to test for numerical unstability                */
/*     x    rwv(ix+1)  : vector containing the abscissae of the starting points  */
/*     hx   rwv(ihx+1) : vector containing the ordinates                         */
/*     hpx  rwv(ihpx+1): vector containing the derivatives                       */
/*     lb   iwv(5) : boolean indicating if there is a lower bound to the domain  */
/*     xlb  rwv(8) : value of the lower bound                                    */
/*     ub   iwv(6) : boolean indicating if there is an upper bound               */
/*     xub  rwv(9) : value of the upper bound                                    */
/*     ifault      : diagnostic                                                  */
/*     iwv,rwv     : integer and real working vectors                            */
/*                                                                               */
/* MORE DETAILED DESCRIPTION OF PARAMETERS                                       */
/*                                                                               */
// parameter  type               i/o     description
// ---------  ----               ---     -----------
//
// SUBROUTINE INITIAL_ (ns, m, emax, x, hx, hpx, lb, xlb, ub, xub, ifault, iwv, rwv)
//
// ns         integer            input   upper limit for number of
//                                       points defining hulls (defines
//                                       lengths of working vectors: see
//                                       below; ns=10 should be ample).     // ample = hojny (CZ)
//
// m          integer            input   number of starting abscissae
//                                       (m=2 is usually enough)
//
// emax       double precision   input   a large number for which
//                                       exp(emax) and exp(-emax) can be
//                                       calculated
//
// x          double precision   input   starting abscissae
//            array size m
//
// hx         double precision   input   log density at starting
//            array size m               abscissae
//
// hpx        double precision   input   gradients of log density at
//            array size m               starting abscissae
//
// lb         logical            input   true if domain is bounded below
//
// xlb        double precision   input   lower bound of domain
//                                      (if lb = true)
//
// ub         logical            input   true if domain is bounded above
//
// xub        double precision   input   upper bound of domain
//                                      (if ub = true)
//
// ifault     integer            output  0: successful initialisation
//                                       1: not enough starting points
//                                       2: ns is less than m
//                                       3: no abscissae to left of mode
//                                          (if lb = false)
//                                       4: no abscissae to right of
//                                          mode (if ub = false)
//                                       5: non-log-concavity detected
//
// iwv        integer array      output  working vector
//            size ns+7
//
// rwv        double precision   output  working vector
//            array size 6*ns+15 = 6*(ns+1)+9
//
void initial_(const int *ns, const int *m, const double *emax, const double *x, const double *hx, const double *hpx,
              const int *lb, double *xlb, const int *ub, double *xub, int *ifault, int *iwv, double *rwv)
{
  bool test = false;

  /* System generated locals */
  static int i__1;
  static double d__1, d__2;

  /* Local variables */
  static double alcu, hulb, huub;
  static int iipt, ihpx, ilow, ihuz, i__, ihigh, iscum;
  static bool horiz;
  static double cu;
  static int nn, ix, iz;
  static double huzmax, eps;
  static int ihx;

  /* Parameter adjustments */
  --rwv;
  --iwv;
  --hpx;
  --hx;
  --x;

  /* Function Body */
  d__1 = -(*emax);
  eps = expon_(&d__1, emax);
  *ifault = 0;
  ilow = 1;
  ihigh = 1;
  nn = *ns + 1;
  /* at least one starting point */
  if (*m < 1)
  {
    *ifault = 1;
  }
  huzmax = hx[1];
  if (!(*ub))
  {
    *xub = (double)0.0;
  }
  if (!(*lb))
  {
    *xlb = (double)0.0;
  }
  hulb = (*xlb - x[1]) * hpx[1] + hx[1];
  huub = (*xub - x[1]) * hpx[1] + hx[1];
  /* if bounded on both sides */
  if (*ub && *lb)
  {
    huzmax = (huub > hulb ? huub : hulb);
    horiz = fabs(hpx[1]) < eps;
    if (horiz)
    {
      d__1 = (huub + hulb) * (double)0.5 - huzmax;
      cu = expon_(&d__1, emax) * (*xub - *xlb);
    }
    else
    {
      d__1 = huub - huzmax;
      d__2 = hulb - huub;
      cu = expon_(&d__1, emax) * (1 - expon_(&d__2, emax)) / hpx[1];
    }
  }
  else if (*ub && !(*lb))
  {
    /* if bounded on the right and unbounded on the left */
    huzmax = huub;
    cu = (double)1.0 / hpx[1];
  }
  else if (!(*ub) && *lb)
  {
    /* if bounded on the left and unbounded on the right */
    huzmax = hulb;
    cu = (double)-1.0 / hpx[1];
    /* if unbounded at least 2 starting points */
  }
  else
  {
    cu = (double)0.0;
    if (*m < 2)
    {
      *ifault = 1;
    }
  }
  if (cu > (double)0.0)
  {
    alcu = log(cu);
  }
  /* set pointers */
  iipt = 6;
  iz = 9;
  ihuz = nn + iz;
  iscum = nn + ihuz;
  ix = nn + iscum;
  ihx = nn + ix;
  ihpx = nn + ihx;
  /* store values in working vectors */
  iwv[1] = ilow;
  iwv[2] = ihigh;
  iwv[3] = *ns;
  iwv[4] = 1;
  if (*lb)
  {
    iwv[5] = 1;
  }
  else
  {
    iwv[5] = 0;
  }
  if (*ub)
  {
    iwv[6] = 1;
  }
  else
  {
    iwv[6] = 0;
  }
  if (*ns < *m)
  {
    *ifault = 2;
  }
  iwv[iipt + 1] = 0;
  rwv[1] = hulb;
  rwv[2] = huub;
  rwv[3] = *emax;
  rwv[4] = eps;
  rwv[5] = cu;
  rwv[6] = alcu;
  rwv[7] = huzmax;
  rwv[8] = *xlb;
  rwv[9] = *xub;
  rwv[iscum + 1] = (double)1.0;
  i__1 = *m;
  for (i__ = 1; i__ <= i__1; ++i__)
  {
    rwv[ix + i__] = x[i__];
    rwv[ihx + i__] = hx[i__];
    rwv[ihpx + i__] = hpx[i__];
    /* L9: */
  }
  /* create lower and upper hulls */
  i__ = 1;
  if (test)
  {
    Rprintf("\n*m = %d, i__ = %d;  ", *m, i__);
    Rprintf("x = ");
    for (int i = 1; i <= 1 + (*ns); i++)
      Rprintf("%f,  ", rwv[ix + i]);
    Rprintf("\n");
    Rprintf("iwv = ");
    for (int i = 1; i <= 6; i++)
      Rprintf("%d,  ", iwv[i]);
    Rprintf("\n");
    Rprintf("ipt = ");
    for (int i = 7; i <= 7 + (*ns); i++)
      Rprintf("%d,  ", iwv[i]);
    Rprintf("\n");
  }
  while (i__ < *m)
  {
    update_(&iwv[4], &iwv[1], &iwv[2], &iwv[iipt + 1], &rwv[iscum + 1],
            &rwv[5], &rwv[ix + 1], &rwv[ihx + 1], &rwv[ihpx + 1], &rwv[iz + 1],
            &rwv[ihuz + 1], &rwv[7], &rwv[3], lb, &rwv[8],
            &rwv[1], ub, &rwv[9], &rwv[2], ifault,
            &rwv[4], &rwv[6]);
    i__ = iwv[4];
    if (test)
    {
      Rprintf("\ni__=%d;  ", i__);
      Rprintf("x = ");
      for (int i = 1; i <= 1 + (*ns); i++)
        Rprintf("%f,  ", rwv[ix + i]);
      Rprintf("\n");
      Rprintf("iwv = ");
      for (int i = 1; i <= 6; i++)
        Rprintf("%d,  ", iwv[i]);
      Rprintf("\n");
      Rprintf("ipt = ");
      for (int i = 7; i <= 7 + (*ns); i++)
        Rprintf("%d,  ", iwv[i]);
      Rprintf("\n");
    }
    if (*ifault != 0)
      return;
  }
  /* test for wrong starting points */
  if (!(*lb) && hpx[iwv[1]] < eps)
  {
    *ifault = 3;
  }
  if (!(*ub) && hpx[iwv[2]] > -eps)
  {
    *ifault = 4;
  }
  return;
} /* end of the routine initial_ */

/* *********************************************************************** */

// Subroutine splhull_:
// ====================
/* This subroutine samples beta from the normalised upper hull */
//
void splhull_(double *u2, int *ipt, int *ilow, const int *lb, double *xlb,
              double *hulb, double *huzmax, double *alcu, double *x, double *hx,
              double *hpx, double *z__, double *huz, double *scum, const double *eps,
              const double *emax, double *beta, int *i__, int *j)
{
  /* System generated locals */
  double d__1, d__2;

  /* Builtin functions */
  //    double log(doublereal);

  /* Local variables */
  static double sign, logdu, logtg;
  static bool horiz;
  static double eh;
  //    extern doublereal expon_(doublereal *, doublereal *);

  /* Parameter adjustments */
  --scum;
  --huz;
  --z__;
  --hpx;
  --hx;
  --x;
  --ipt;

  /* Function Body */
  *i__ = *ilow;

  /* find from which exponential piece you sample */
  while (*u2 > scum[*i__])
  {
    *j = *i__;
    *i__ = ipt[*i__];
  }
  if (*i__ == *ilow)
  {
    /* sample below z(ilow),depending on the existence of a lower bound */
    if (*lb)
    {
      eh = *hulb - *huzmax - *alcu;
      horiz = (d__1 = hpx[*ilow], fabs(d__1)) < *eps;
      if (horiz)
      {
        d__1 = -eh;
        *beta = *xlb + *u2 * expon_(&d__1, emax);
      }
      else
      {
        sign = (d__1 = hpx[*i__], fabs(d__1)) / hpx[*i__];
        d__2 = (d__1 = hpx[*i__], fabs(d__1));
        logtg = log(d__2);
        logdu = log(*u2);
        eh = logdu + logtg - eh;
        if (eh < *emax)
        {
          d__1 = sign * expon_(&eh, emax) + (double)1.0;
          *beta = *xlb + log(d__1) / hpx[*i__];
        }
        else
        {
          *beta = *xlb + eh / hpx[*i__];
        }
      }
    }
    else
    {
      /*     hpx(i) must be positive , x(ilow) is left of the mode */
      d__1 = hpx[*i__] * *u2;
      *beta = (log(d__1) + *alcu - hx[*i__] + x[*i__] * hpx[*i__] + *huzmax) / hpx[*i__];
    }
  }
  else
  {
    /*   sample above(j) */
    eh = huz[*j] - *huzmax - *alcu;
    horiz = (d__1 = hpx[*i__], fabs(d__1)) < *eps;
    if (horiz)
    {
      d__1 = -eh;
      *beta = z__[*j] + (*u2 - scum[*j]) * expon_(&d__1, emax);
    }
    else
    {
      sign = (d__1 = hpx[*i__], fabs(d__1)) / hpx[*i__];
      d__2 = (d__1 = hpx[*i__], fabs(d__1));
      logtg = log(d__2);
      d__1 = *u2 - scum[*j];
      logdu = log(d__1);
      eh = logdu + logtg - eh;
      if (eh < *emax)
      {
        d__1 = sign * expon_(&eh, emax) + (double)1.0;
        *beta = z__[*j] + log(d__1) / hpx[*i__];
      }
      else
      {
        *beta = z__[*j] + eh / hpx[*i__];
      }
    }
  }

  return;
} /* end of the routine splhull_ */

/* *********************************************************************** */

// Subroutine intersection_:
// =========================
/* computes the intersection (z1,hz1) between 2 tangents defined by */
/*   x1,y1,yp1 and x2,y2,yp2 */
//
void intersection_(double *x1, double *y1, double *yp1, double *x2, double *y2,
                   double *yp2, double *z1, double *hz1, const double *eps, int *ifault)
{
  static double dh, y12, y21;

  /* first test for non-concavity */
  y12 = *y1 + *yp1 * (*x2 - *x1);
  y21 = *y2 + *yp2 * (*x1 - *x2);
  if (y21 < *y1 || y12 < *y2)
  {
    //      REprintf("\nTrap: non-logcocavity detected by ARS intersection_ function\ny21=%15.15e, y12=%15.15e\n", y21, y12);
    //      REprintf("*x1=%15.10e, *x2=%15.10e, *y1=%15.15e, *y2=%15.15e, *yp1=%15.10e, *yp2=%15.10e\n", *x1, *x2, *y1, *y2, *yp1, *yp2);
    //      if (y21 < *y1) REprintf("y21 < *y1\n");
    //      if (y12 < *y2) REprintf("y12 < *y2\n");
    *ifault = 5;
    return;
  }
  dh = *yp2 - *yp1;
  /*  IF the lines are nearly parallel, */
  /*  the intersection is taken at the midpoint */
  if (fabs(dh) <= *eps)
  {
    *z1 = (*x1 + *x2) * (double)0.5;
    *hz1 = (*y1 + *y2) * (double)0.5;
    /*  Else compute from the left or the right for greater numerical */
    /*       precision */
  }
  else if (fabs(*yp1) < fabs(*yp2))
  {
    *z1 = *x2 + (*y1 - *y2 + *yp1 * (*x2 - *x1)) / dh;
    *hz1 = *yp1 * (*z1 - *x1) + *y1;
  }
  else
  {
    *z1 = *x1 + (*y1 - *y2 + *yp2 * (*x2 - *x1)) / dh;
    *hz1 = *yp2 * (*z1 - *x2) + *y2;
  }
  /*  test for misbehaviour due to numerical imprecision */
  if (*z1 < *x1 || *z1 > *x2)
  {
    *ifault = 7;
  }
  return;
} /* end of the routine intersection_ */

/* *********************************************************************** */

// Subroutine update_:
// ===================
/* This subroutine increments n and updates all the parameters which */
/* define the lower and the upper hull */
//
/* DESCRIPTION OF PARAMETERS and place of storage */

/*     ilow iwv(1)    : index of the smallest x(i) */
/*     ihigh iwv(2)   : index of the largest x(i) */
/*     n    iwv(4)    : number of points defining the hulls */
/*     ipt  iwv(iipt) : pointer array:  ipt(i) is the index of the x(.) */
/*                      immediately larger than x(i) */
/*     hulb rwv(1)    : value of the upper hull at xlb */
/*     huub rwv(2)    : value of the upper hull at xub */
/*     cu   rwv(5)    : integral of the exponentiated upper hull divided */
/*                      by exp(huzmax) */
/*     alcu rwv(6)    : logarithm of cu */
/*     huzmax rwv(7)  : maximum of huz(i); i=1,n */
/*     z    rwv(iz+1) : z(i) is the abscissa of the intersection between */
/*                      the tangents at x(i) and x(ipt(i)) */
/*     huz  rwv(ihuz+1): huz(i) is the ordinate of the intersection */
/*                        defined above */
/*     scum rwv(iscum): scum(i) is the cumulative probability of the */
/*                      normalised exponential of the upper hull */
/*                      calculated at z(i) */
/*     eps  rwv(4)    : =exp(-emax) a very small number */
//
void update_(int *n, int *ilow, int *ihigh, int *ipt, double *scum,
             double *cu, double *x, double *hx, double *hpx, double *z__,
             double *huz, double *huzmax, const double *emax,
             const int *lb, double *xlb, double *hulb, const int *ub, double *xub, double *huub,
             int *ifault, const double *eps, double *alcu)
{
  const double ZERO_2_DERIV = 1e-2; // zero for test performed by this routine
                                    //  * this should detect non-zero second derivative
                                    //  * higher values are more safe and we only loose very small amount of time
  //    bool test = false;

  /* System generated locals */
  double d__1, d__2;

  /* Local variables */
  static int i__, j;
  static double u;
  static bool horiz;
  static double dh;

  /* Parameter adjustments */
  --huz;
  --z__;
  --hpx;
  --hx;
  --x;
  --scum;
  --ipt;

  /* Function Body */
  ++(*n);
  /* update z,huz and ipt */
  if (x[*n] < x[*ilow])
  {
    /* insert x(n) below x(ilow) */
    /*   test for non-concavity */
    if (hpx[*ilow] > hpx[*n])
    {
      //	  REprintf("Trap: non-logcocavity detected by ARS update_ function\nhpx[*ilow]=%e, hpx[*n]=%e\n", hpx[*ilow], hpx[*n]);
      *ifault = 5;
    }
    ipt[*n] = *ilow;
    intersection_(&x[*n], &hx[*n], &hpx[*n], &x[*ilow], &hx[*ilow],
                  &hpx[*ilow], &z__[*n], &huz[*n], eps, ifault);
    if (*ifault != 0)
    {
      return;
    }
    if (*lb)
    {
      *hulb = hpx[*n] * (*xlb - x[*n]) + hx[*n];
    }
    *ilow = *n;
  }
  else
  {
    i__ = *ilow;
    j = i__;
    /* find where to insert x(n) */
    while (x[*n] >= x[i__] && ipt[i__] != 0)
    {
      j = i__;
      i__ = ipt[i__];
    }
    if (x[*n] >= x[i__])
    {
      /* insert above x(ihigh) */
      /*   test for non-concavity */
      if (hpx[i__] < hpx[*n])
      {
        //   	        REprintf("Trap: non-logcocavity detected by ARS update_ function\nhpx[i__]=%e, hpx[*n]=%e\n", hpx[i__], hpx[*n]);
        *ifault = 5;
      }
      *ihigh = *n;
      ipt[i__] = *n;
      ipt[*n] = 0;
      intersection_(&x[i__], &hx[i__], &hpx[i__], &x[*n], &hx[*n],
                    &hpx[*n], &z__[i__], &huz[i__], eps, ifault);
      if (*ifault != 0)
      {
        return;
      }
      *huub = hpx[*n] * (*xub - x[*n]) + hx[*n];
      z__[*n] = (double)0.0;
      huz[*n] = (double)0.0;
    }
    else
    {
      /* insert x(n) between x(j) and x(i) */
      /*   test for non-concavity */
      if (hpx[j] < hpx[*n] || hpx[i__] > hpx[*n])
      {
        //	      REprintf("Trap: non-logcocavity detected by ARS update_ function\nhpx[j]=%e, hpx[i__]=%e, hpx[*n]=%e\n", hpx[j], hpx[i__], hpx[*n]);
        *ifault = 5;
      }
      ipt[j] = *n;
      ipt[*n] = i__;
      /*     insert z(j) between x(j) and x(n) */
      intersection_(&x[j], &hx[j], &hpx[j], &x[*n], &hx[*n],
                    &hpx[*n], &z__[j], &huz[j], eps, ifault);
      if (*ifault != 0)
      {
        return;
      }
      /*     insert z(n) between x(n) and x(i) */
      intersection_(&x[*n], &hx[*n], &hpx[*n], &x[i__], &hx[i__],
                    &hpx[i__], &z__[*n], &huz[*n], eps, ifault);
      if (*ifault != 0)
      {
        return;
      }
    }
  }
  /* update huzmax */
  j = *ilow;
  i__ = ipt[j];
  *huzmax = huz[j];
  while (huz[j] < huz[i__] && ipt[i__] != 0)
  {
    j = i__;
    i__ = ipt[i__];
    /* Computing MAX */
    d__1 = *huzmax, d__2 = huz[j];
    *huzmax = (d__1 > d__2 ? d__1 : d__2);
  }
  if (*lb)
  {
    *huzmax = (*huzmax > *hulb ? *huzmax : *hulb);
  }
  if (*ub)
  {
    *huzmax = (*huzmax > *huub ? *huzmax : *huub);
  }
  /* update cu */
  /*  scum receives area below exponentiated upper hull left of z(i) */
  i__ = *ilow;
  horiz = (d__1 = hpx[*ilow], fabs(d__1)) < *eps;
  if (!(*lb) && !horiz)
  {
    d__1 = huz[i__] - *huzmax;
    *cu = expon_(&d__1, emax) / hpx[i__];
  }
  else if (*lb && horiz)
  {
    d__1 = *hulb - *huzmax;
    *cu = (z__[*ilow] - *xlb) * expon_(&d__1, emax);
  }
  else if (*lb && !horiz)
  {
    dh = *hulb - huz[i__];
    if (dh > *emax)
    {
      d__1 = *hulb - *huzmax;
      *cu = -expon_(&d__1, emax) / hpx[i__];
    }
    else
    {
      d__1 = huz[i__] - *huzmax;
      *cu = expon_(&d__1, emax) * (1 - expon_(&dh, emax)) / hpx[i__];
    }
  }
  else
  {
    *cu = 0.;
  }
  scum[i__] = *cu;
  j = i__;
  i__ = ipt[i__];

  int control_count = 0; /* added by AK */
  while (ipt[i__] != 0)
  {
    if (control_count > *n)
      stop("Trap in ARS: infinite while in update_ of ars.cpp near l. 810\n"); /* added by AK */
    control_count++;                                                           /* added by AK */
    dh = huz[j] - huz[i__];
    horiz = (d__1 = hpx[i__], fabs(d__1)) < *eps;
    if (horiz)
    {
      d__1 = (huz[i__] + huz[j]) * (double)0.5 - *huzmax;
      *cu += (z__[i__] - z__[j]) * expon_(&d__1, emax);
    }
    else
    {
      if (dh < *emax)
      {
        d__1 = huz[i__] - *huzmax;
        *cu += expon_(&d__1, emax) * (1 - expon_(&dh, emax)) / hpx[i__];
      }
      else
      {
        d__1 = huz[j] - *huzmax;
        *cu -= expon_(&d__1, emax) / hpx[i__];
      }
    }
    j = i__;
    i__ = ipt[i__];
    scum[j] = *cu;
  }
  horiz = (d__1 = hpx[i__], fabs(d__1)) < *eps;
  /* if the derivative is very small the tangent is nearly horizontal */
  if (!(*ub || horiz))
  {
    d__1 = huz[j] - *huzmax;
    *cu -= expon_(&d__1, emax) / hpx[i__];
  }
  else if (*ub && horiz)
  {
    d__1 = (*huub + hx[i__]) * (double)0.5 - *huzmax;
    *cu += (*xub - x[i__]) * expon_(&d__1, emax);
  }
  else if (*ub && !horiz)
  {
    dh = huz[j] - *huub;
    if (dh > *emax)
    {
      d__1 = huz[j] - *huzmax;
      *cu -= expon_(&d__1, emax) / hpx[i__];
    }
    else
    {
      d__1 = *huub - *huzmax;
      *cu += expon_(&d__1, emax) * (1 - expon_(&dh, emax)) / hpx[i__];
    }
  }
  scum[i__] = *cu;
  if (*cu > 0.)
  {
    *alcu = log(*cu);
  }
  /* normalize scum to obtain a cumulative probability while excluding */
  /*    unnecessary points */
  i__ = *ilow;
  u = (*cu - scum[i__]) / *cu;
  if (u == (double)1.0 && hpx[ipt[i__]] > ZERO_2_DERIV)
  {
    *ilow = ipt[i__];
    scum[i__] = (double)0.0;
  }
  else
  {
    scum[i__] = (double)1.0 - u;
  }
  j = i__;
  i__ = ipt[i__];
  while (ipt[i__] != 0)
  {
    j = i__;
    i__ = ipt[i__];
    u = (*cu - scum[j]) / *cu;
    if (u == (double)1.0 && hpx[i__] > ZERO_2_DERIV)
    {
      *ilow = i__;
    }
    else
    {
      scum[j] = (double)1.0 - u;
    }
  }
  scum[i__] = (double)1.0;
  if (*ub)
  {
    *huub = hpx[*ihigh] * (*xub - x[*ihigh]) + hx[*ihigh];
  }
  if (*lb)
  {
    *hulb = hpx[*ilow] * (*xlb - x[*ilow]) + hx[*ilow];
  }
  return;
} /* end of the routine update_ */
