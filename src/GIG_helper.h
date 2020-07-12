#ifndef GIG_HELPER_H
#define GIG_HELPER_H

double _gig_mode(double lambda, double omega);
void _rgig_ROU_noshift (double *res, int n, double lambda, double lambda_old, double omega, double alpha);
void _rgig_newapproach1 (double *res, int n, double lambda, double lambda_old, double omega, double alpha);
void _rgig_ROU_shift_alt (double *res, int n, double lambda, double lambda_old, double omega, double alpha);
double rgig(double lambda, double chi, double psi);
  



#endif