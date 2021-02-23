// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
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
#include "ars_helper.h"
//#include "ars_pois_helper.h"

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

void spl1_(const int *ns, int *n, int *ilow, int *ihigh, int *ipt,
           double *scum, double *cu, double *x, double *hx, double *hpx,
           double *z__, double *huz, double *huzmax,
           const int *lb, double *xlb, double *hulb, const int *ub, double *xub, double *huub,
           int *ifault, const double *emax, const double *eps, double *alcu,
           int l, int w, // which node
           arma::mat &Z_curr,
           const arma::mat &mu_Z,
           const arma::mat &Sigma_Z, // this is Sigma (cov) not Omega (percision)
           const arma::mat &y,
           int k, int p, int n_sample)
{
  const int max_attempt = 10 * (*ns); // maximal number of attempts to sample a value
  //  (usually (not necessarily) is something wrong if this number is reached)

  /* Local variables */
  static double alhl, alhu;
  static int i__, j, n1;
  static double u1, u2, fx;
  static bool sampld;
  static double alu1;

  double temp = 0; // to save new sample of Z(l,w)
  double *beta = &temp;

  //Necesario para poder utilizar los generadores de numeros aleatorios del R
  GetRNGstate();

  /* Parameter adjustments */
  --huz;
  --z__;
  --hpx;
  --hx;
  --x;
  --scum;
  --ipt;

  // calculate posterior related things
  //arma::vec res(2,fill::zeros);
  //hx += (y(i,j) * Z_curr(i,j) - exp(Z_curr(i,j))); // log posterior due to Poisson
  //res(1) += (y(i,j) - exp(Z_curr(i,j))); // d/dz logPost due to Poisson
  //Rcout << res <<endl;
  arma::uvec ind = linspace<uvec>(0, k - 1, k);
  arma::uvec indi = linspace<uvec>(0, n_sample - 1, n_sample);
  arma::uvec ind_noj = find(ind != w);

  double Sigmabb = Sigma_Z(w, w);
  arma::mat Sigmac = Sigma_Z(ind_noj, find(ind == w));
  arma::mat Sigmaa = Sigma_Z(ind_noj, ind_noj);

  double mu_Zij = mu_Z(l, w) + as_scalar(trans(Sigmac) * solve(Sigmaa, trans(Z_curr(find(indi == l), ind_noj) - mu_Z(find(indi == l), ind_noj))));
  double sigma2_Zij = Sigmabb - as_scalar(trans(Sigmac) * solve(Sigmaa, Sigmac));
  //mu_Zij = as_scalar(mu_Zij);
  //sigma2_Zij = as_scalar(sigma_Zij);

  //double h = 0 ;//= (0.5 * (*beta-mu_Zij)*(Z_curr(i,j)-mu_Zij)/sigma2_Zij)+(y(i,j) * (*beta) - exp(*beta));
  //double hprime = 0; //  = ((*beta-mu_Zij)/sigma2_Zij)+(y(i,j) - exp(*beta));

  /* Function Body */
  *ifault = 0;
  sampld = false;
  int attempts = 0;
  while (!sampld && attempts < max_attempt)
  {
    //	u2 = random_(&l);
    u2 = unif_rand();
    /* test for zero random number */
    if (u2 == (double)0.0)
    {
      *ifault = 6;
      return;
    }
    splhull_(&u2, &ipt[1], ilow, lb, xlb,
             hulb, huzmax, alcu, &x[1], &hx[1],
             &hpx[1], &z__[1], &huz[1], &scum[1], eps,
             emax, beta, &i__, &j);
    /* sample u1 to compute rejection */
    u1 = unif_rand();
    if (u1 == (double)0.0)
    {
      *ifault = 6;
    }
    alu1 = log(u1);
    /* compute alhu: upper hull at point u1 */
    alhu = hpx[i__] * (*beta - x[i__]) + hx[i__] - *huzmax;
    if (*beta > x[*ilow] && *beta < x[*ihigh])
    {
      /* compute alhl: value of the lower hull at point u1 */
      if (*beta > x[i__])
      {
        j = i__;
        i__ = ipt[i__];
      }
      alhl = hx[i__] + (*beta - x[i__]) * (hx[i__] - hx[j]) / (x[i__] - x[j]) - *huzmax;
      /* squeezing test */
      if (alhl - alhu > alu1)
      {
        sampld = true;
      }
      //            else{
      //              Rprintf("alhl=%e, alhu=%e, alu1=%e\n", alhl, alhu, alu1);
      //            }
    }
    /* if not sampled evaluate the function, do the rejection test and update */
    if (!sampld)
    {
      n1 = *n + 1;
      x[n1] = *beta;
      // h and hprime
      hx[n1] = (-0.5 * (*beta - mu_Zij) * (*beta - mu_Zij) / sigma2_Zij) + (y(l, w) * (*beta) - exp(*beta));

      hpx[n1] = (-(*beta - mu_Zij) / sigma2_Zij) + (y(l, w) - exp(*beta));
      ;
      fx = hx[n1] - *huzmax;
      if (alu1 < fx - alhu)
      {
        sampld = true;
      }
      //            else{
      //              Rprintf("alu1=%e, fx=%e, alhu=%e\n", alu1, fx, alhu);
      //            }
      /* update while the number of points defining the hulls is lower than ns */
      if (*n < *ns)
      {
        update_(n, ilow, ihigh, &ipt[1], &scum[1],
                cu, &x[1], &hx[1], &hpx[1], &z__[1],
                &huz[1], huzmax, emax, lb, xlb,
                hulb, ub, xub, huub, ifault,
                eps, alcu);
      }
      if (*ifault != 0)
      {
        return;
      }
    }
    attempts++;
  } /** end of while (! sampld) **/
  //Necesario al terminar de utilizar los generadores de numeros aleatorios del R
  PutRNGstate();
  //if (attempts >= max_attempt)
  //  Rcout << "Trap in ARS: Maximum number of attempts reached by routine spl1_\n"
  //        << endl;
  Z_curr(l, w) = *beta;
  return;
} /* end of the routine spl1_ */

void sample_(int *iwv, double *rwv,
             int i, int j, // which node
             arma::mat &Z_curr,
             const arma::mat &mu_Z,
             const arma::mat &Sigma_Z, // this is Sigma (cov) not Omega (percision)
             const arma::mat &y,
             int k, int p, int n, int *ifault)
{
  static int iipt, ihpx, ihuz, iscum;
  static int lb, ub;
  static int nn, ns, ix, iz, ihx;

  /* Parameter adjustments */
  --rwv;
  --iwv;

  /* Function Body */
  iipt = 6;
  iz = 9;
  ns = iwv[3];
  nn = ns + 1;
  ihuz = nn + iz;
  iscum = nn + ihuz;
  ix = nn + iscum;
  ihx = nn + ix;
  ihpx = nn + ihx;
  lb = 0; // false
  ub = 0; // false
  if (iwv[5] == 1)
  {
    lb = 1; // true
  }
  if (iwv[6] == 1)
  {
    ub = 1; // true
  }

  /*     call sampling subroutine */
  spl1_(&ns, &iwv[4], &iwv[1], &iwv[2], &iwv[iipt + 1],
        &rwv[iscum + 1], &rwv[5], &rwv[ix + 1], &rwv[ihx + 1], &rwv[ihpx + 1],
        &rwv[iz + 1], &rwv[ihuz + 1], &rwv[7], &lb, &rwv[8],
        &rwv[1], &ub, &rwv[9], &rwv[2], ifault,
        &rwv[3], &rwv[4], &rwv[6], i, j, // which node
        Z_curr,
        mu_Z,
        Sigma_Z, // this is Sigma (cov) not Omega (percision)
        y,
        k, p, n);
  return;
} /* end of the routine sample_ */

// [[Rcpp::export]]
void update_Z_helper_Pois(arma::mat &Z_curr,
                          const arma::mat &mu_Z,
                          const arma::mat &Sigma_Z, // this is Sigma (cov) not Omega (percision)
                          const arma::mat &y,
                          int k, int p, int n,
                          int ns, int m, double emax // ars parameters
)
{

  arma::uvec ind = linspace<uvec>(0, k - 1, k);
  arma::uvec indi = linspace<uvec>(0, n - 1, n);
  arma::uvec ind_noj;
  double Sigmabb, mu_Zij, sigma2_Zij;
  arma::mat Sigmac, Sigmaa;

  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < k; ++j)
    {

      ind_noj = find(ind != j);

      Sigmabb = Sigma_Z(j, j);
      Sigmac = Sigma_Z(ind_noj, find(ind == j));
      Sigmaa = Sigma_Z(ind_noj, ind_noj);

      mu_Zij = mu_Z(i, j) + as_scalar(trans(Sigmac) * solve(Sigmaa, trans(Z_curr(find(indi == i), ind_noj) - mu_Z(find(indi == i), ind_noj))));
      sigma2_Zij = Sigmabb - as_scalar(trans(Sigmac) * solve(Sigmaa, Sigmac));

      int *iwv = new int[ns + 7]();
      double *rwv = new double[6 * (ns + 1) + 9]();
      //double *x = new double[ns]();
      //double *hx = new double[ns]();
      //double *hpx = new double[ns]();

      double *x = new double[m]();
      double *hx = new double[m]();
      double *hpx = new double[m]();
      vec h_hprime_temp(2);

      int lb = 0;
      int ub = 0;
      double xlb = 0;
      double xub = 0;
      int ifault = 0;

      double center = (log(y(i, j) + .01) + mu_Zij) / 2 ;
      double range = sqrt(sigma2_Zij);
      bool bad_init = true;
      double left_hp, right_hp;
      while(bad_init){
         left_hp = -((center - range - mu_Zij) / sigma2_Zij) + (y(i, j) - exp(center-range));
         right_hp = -((center + range - mu_Zij) / sigma2_Zij) + (y(i, j) - exp(center+range));
         bad_init = left_hp * right_hp >= 0;
         if(bad_init) {range += sqrt(sigma2_Zij);}// adaptively chose intial points
      }
      range *= 1.1;// being safe


      //Rcout<< "before ars" << i << " " << j << "\n" << Z_curr(i,j) <<endl;
      for (int ww = 0; ww < m; ++ww)
      {
        x[ww] = center + ((double)ww - ((double)m / 2)) * (2*range/m);
        //Rcout << log(y(i,j)+.01) + ((double)ww-((double)m/2)) * (4/(double)m) << "  " << x[ww] <<endl;
        //Z_curr(i,j) = x[ww];
        //Rcout << "ars working" <<endl;
        //Rcout << "i:" << i << "  j:" << j << "  ww:" << ww <<endl;
        //Rcout << "flag" <<endl;

        //Rcout << h_hprime_temp <<endl;
        hx[ww] = -(0.5 * (x[ww] - mu_Zij) * (x[ww] - mu_Zij) / sigma2_Zij) + (y(i, j) * (x[ww]) - exp(x[ww]));

        hpx[ww] = -((x[ww] - mu_Zij) / sigma2_Zij) + (y(i, j) - exp(x[ww]));

        //cout << x[ww] << "  " << hx[ww] << "  " << hpx[ww] <<endl;
      } //initial support

      initial_(&ns, &m, &emax, x, hx, hpx,
               &lb, &xlb, &ub, &xub, &ifault, iwv, rwv);
      //for(int kkk = 0 ; kkk < 7 ; kkk++) Rcout << "iwv "<< kkk <<": " << iwv[kkk]<<" ";
      //Rcout<<endl;
      //Rcout<< "i: " << i << " j: " << j << " y:  " << y(i,j)<< "  sigma2: " << sigma2_Zij << " mu: " << mu_Zij << endl;
      sample_(iwv, rwv,
              i, j, // which node
              Z_curr,
              mu_Z,
              Sigma_Z, // this is Sigma (cov) not Omega (percision)
              y,
              k, p, n, &ifault);
      //Rcout << "after:\n" << Z_curr(i,j) <<endl;
      //if (ifault != 0)
      //{
        //Rcout << "ARS failed with code" << ifault << endl;
        //stop("ARS failed with code %i \n",ifault);
      //}
      delete[] iwv;
      delete[] rwv;
      delete[] x;
      delete[] hx;
      delete[] hpx;
    }
  }

  //Rcout << Z_curr <<endl;
  return;
}

// [[Rcpp::export]]
void update_Z_helper_Pois_reg(arma::mat &Z_curr, // persumably large, thus will not copy
                              const arma::mat &data,
                              const arma::mat &design,
                              const arma::vec &mu_curr,
                              const arma::mat &beta_curr,
                              const arma::mat &Omega_curr,
                              int k, int p, int n,
                              int ns, int m, double emax // ars parameters
)
{
  arma::mat mu_Zmat = design * beta_curr;
  mu_Zmat.each_row() += mu_curr.t(); // calculate the expectation of latent
  arma::mat Sigma_Z = inv_sympd(Omega_curr);
  //Rcout << "muZ:\n" << mu_Zmat <<endl;
  //Rcout << "SigmaZ:\n" << Sigma_Z <<endl;
  update_Z_helper_Pois(Z_curr, mu_Zmat, Sigma_Z, data,
                       k, p, n, ns, m, emax);

  //Rcout <<Z_curr<<endl;
  return;
}

void update_Z_helper_Pois_gra(arma::mat &Z_curr, // persumably large, thus will not copy
                              const arma::mat &data,
                              const arma::vec &mu_curr,
                              const arma::mat &Omega_curr,
                              int k, int p, int n,
                              int ns, int m, double emax // ars parameters
)
{
  arma::mat mu_Zmat = 0 * Z_curr;
  mu_Zmat.each_row() += mu_curr.t(); // calculate the expectation of latent
  arma::mat Sigma_Z = inv_sympd(Omega_curr);
  update_Z_helper_Pois(Z_curr, mu_Zmat, Sigma_Z, data,
                       k, p, n, ns, m, emax);
  return;
}

// [[Rcpp::export]]
void update_Z_helper_Pois_CAR(arma::mat &Z_curr, // persumably large, thus will not copy
                              const arma::mat &data,
                              const arma::mat &design,
                              const arma::vec &mu_curr,
                              const arma::mat &beta_curr,
                              const arma::mat &Omega_curr,
                              int k, int p, int n,
                              int ns, int m, double emax // ars parameters
)
{
  arma::mat mu_Zmat = design * beta_curr;
  mu_Zmat.each_row() += mu_curr.t(); // calculate the expectation of latent
  arma::mat Sigma_Z = inv_sympd(Omega_curr);
  mu_Zmat = mu_Zmat * Sigma_Z; // slightly different from regression, CAR need to times Sigma to mu
  update_Z_helper_Pois(Z_curr, mu_Zmat, Sigma_Z, data,
                       k, p, n, ns, m, emax);
  return;
}

// try parallel version of update_Z_helper_multinomial()


// [[Rcpp::export]]
void update_Z_helper_Pois_CAR_randeff(arma::mat &Z_curr, // persumably large, thus will not copy
                                     const arma::mat &data,
                                     const arma::mat &design,
                                     const arma::mat &design_r,
                                     const arma::vec &mu_curr,
                                     const arma::mat &beta_curr,
                                     const arma::mat &nu_curr,
                                     const arma::mat &Omega_curr,
                                     int k, int p, int n,
                                     int ns, int m, double emax // ars parameters
)
{
  arma::mat mu_Zmat = design * beta_curr + design_r * nu_curr;
  mu_Zmat.each_row() += mu_curr.t(); // calculate the expectation of latent
  arma::mat Sigma_Z = inv_sympd(Omega_curr);
  mu_Zmat = mu_Zmat * Sigma_Z; // slightly different from regression, CAR need to times Sigma to mu
  update_Z_helper_Pois(Z_curr, mu_Zmat, Sigma_Z, data,
                              k, p, n, ns, m, emax);
  return;
}



// struct get_Z_worker : public Worker
// {
//   arma::mat &Z_curr;
//   const arma::mat &mu_Z;
//   const arma::mat &Sigma_Z; // this is Sigma (cov) not Omega (percision)
//   const arma::mat &y;
//   const int &k;
//   const int &p;
//   const int &n;
//   const int &ns;
//   const int &m;
//   const double &emax; // ars parameters

//   get_Z_worker(arma::mat &Z_curr,
//                const arma::mat &mu_Z,
//                const arma::mat &Sigma_Z, // this is Sigma (cov) not Omega (percision)
//                const arma::mat &y,
//                const int &k,
//                const int &p,
//                const int &n,
//                const int &ns,
//                const int &m,
//                const double &emax) : Z_curr(Z_curr), mu_Z(mu_Z), Sigma_Z(Sigma_Z),
//                                      y(y), k(k), p(p), n(n), ns(ns), m(m), emax(emax)
//   {
//   }

//   void operator()(std::size_t begin, std::size_t end)
//   {
//     arma::uvec ind = linspace<uvec>(0, k - 1, k);
//     arma::uvec indi = linspace<uvec>(0, n - 1, n);
//     arma::uvec ind_noj;
//     double Sigmabb, mu_Zij, sigma2_Zij;
//     arma::mat Sigmac, Sigmaa;
//     arma::mat Z__j; // Z w/o the column focuing on
//     double normalizingwoZi;
//     arma::vec N = sum(y, 1);

//     for (int i = begin; i < end; ++i)
//     {
//       for (int j = 0; j < k; ++j)
//       {

//         ind_noj = find(ind != j);

//         Sigmabb = Sigma_Z(j, j);
//         Sigmac = Sigma_Z(ind_noj, find(ind == j));
//         Sigmaa = Sigma_Z(ind_noj, ind_noj);

//         mu_Zij = mu_Z(i, j) + as_scalar(trans(Sigmac) * solve(Sigmaa, trans(Z_curr(find(indi == i), ind_noj) - mu_Z(find(indi == i), ind_noj))));
//         sigma2_Zij = Sigmabb - as_scalar(trans(Sigmac) * solve(Sigmaa, Sigmac));

//         int *iwv = new int[ns + 7]();
//         double *rwv = new double[6 * (ns + 1) + 9]();
//         //double *x = new double[ns]();
//         //double *hx = new double[ns]();
//         //double *hpx = new double[ns]();

//         double *x = new double[m]();
//         double *hx = new double[m]();
//         double *hpx = new double[m]();
//         vec h_hprime_temp(2);

//         int lb = 0;
//         int ub = 0;
//         double xlb = 0;
//         double xub = 0;
//         int ifault = 0;
//         //Rcout<< "before ars" << i << " " << j << "\n" << Z_curr(i,j) <<endl;
//         for (int ww = 0; ww < m; ++ww)
//         {
//           x[ww] = (log(y(i, j) + .01) + mu_Zij) / 2 + ((double)ww - ((double)m / 2)) * ((8 * sqrt(sigma2_Zij)) / (double)m);
//           //Rcout << log(y(i,j)+.01) + ((double)ww-((double)m/2)) * (4/(double)m) << "  " << x[ww] <<endl;
//           //Z_curr(i,j) = x[ww];
//           //Rcout << "ars working" <<endl;
//           //Rcout << "i:" << i << "  j:" << j << "  ww:" << ww <<endl;
//           //Rcout << "flag" <<endl;

//           //Rcout << h_hprime_temp <<endl;
//           hx[ww] = -(0.5 * (x[ww] - mu_Zij) * (x[ww] - mu_Zij) / sigma2_Zij) + (y(i, j) * (x[ww]) - exp(x[ww]));

//           hpx[ww] = -((x[ww] - mu_Zij) / sigma2_Zij) + (y(i, j) - exp(x[ww]));

//           //cout << x[ww] << "  " << hx[ww] << "  " << hpx[ww] <<endl;
//         } //initial support

//         initial_(&ns, &m, &emax, x, hx, hpx,
//                  &lb, &xlb, &ub, &xub, &ifault, iwv, rwv);
//         //for(int kkk = 0 ; kkk < 7 ; kkk++) Rcout << "iwv "<< kkk <<": " << iwv[kkk]<<" ";
//         //Rcout<<endl;
//         //Rcout<< "i: " << i << " j: " << j << " y:  " << y(i,j)<< "  sigma2: " << sigma2_Zij << " mu: " << mu_Zij << endl;
//         sample_(iwv, rwv,
//                 i, j, // which node
//                 Z_curr,
//                 mu_Z,
//                 Sigma_Z, // this is Sigma (cov) not Omega (percision)
//                 y,
//                 k, p, n, &ifault);
//         //Rcout << "after:\n" << Z_curr(i,j) <<endl;
//         if (ifault != 0)
//         {
//           Rcout << "ARS failed with code" << ifault << endl;
//           //stop("ARS failed with code %i \n",ifault);
//         }
//         delete[] iwv;
//         delete[] rwv;
//         delete[] x;
//         delete[] hx;
//         delete[] hpx;
//       }
//     }
//   }
// };

// // [[Rcpp::export]]
// void update_Z_helper_Pois_para(arma::mat &Z_curr,
//                                       const arma::mat &mu_Z,
//                                       const arma::mat &Sigma_Z, // this is Sigma (cov) not Omega (percision)
//                                       const arma::mat &y,
//                                       int k, int p, int n,
//                                       int ns, int m, double emax // ars parameters
// )
// {

//   get_Z_worker Z_worker(Z_curr, mu_Z, Sigma_Z, y,
//                         k, p, n, ns, m, emax);

//   parallelFor(0, n, Z_worker);
//   return;
// }