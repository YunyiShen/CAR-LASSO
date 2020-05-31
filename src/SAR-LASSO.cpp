// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
static const double pi = 3.141592653589793238462643383280;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// Auto-normal

/* Starting here is SAR-LASSO and convert it to CAR if needed, see Ver Hoef et al. 2018
 * The Reason for using SAR rather than CAR is that SAR asks for a ralative easier to check condition 
 *   condition on B (I-B) non-singular, while CAR asks for C to have positive eigen value
 *   We can use independent Laplace prior in SAR to get similar result as LASSO.
 * This is all because Sigma has to be positive semi-definite.
 * Also, matrix B can be asymmertic, 
 *   however this will cause some explaning problem, and identification problem 
 *    as B is no longer an adjacency matrix of the faithfuil graph related with the joint distribution 
 * However, SAR parameterization allows us to do regressive LASSO and graphical LASSO at the same time 
 *   by setting the "noise" in auto regression being iid, this can approximate any graph 
 */

arma::mat getBmat_helper(const arma::vec &B, int k){
  arma::mat Bmat(k,k,fill::zeros);
  arma::vec curr_col(k-1,fill::zeros);
  Bmat(span(1,k-1),0) = B(span(0,k-2)); //first column
  Bmat(span(0,k-2),k-1) = B(span((k-1) * (k-1),k * (k-1)-1)); // last column
  
  for(int i = 1 ; i < k-1 ; ++i){ // columns in between
    curr_col = B( span(i * (k-1), (i+1)*(k-1)-1 ) );
    Bmat(span(0,i-1),i) = curr_col(span(0,i-1));
    Bmat(span(i+1,k-1),i) = curr_col(span(i,k-2));
  }
  return(Bmat);
}



// [[Rcpp::export]]
double Loglik_SAR_Cpp(const arma::vec & Z,
                      const arma::vec & mu,
                      const arma::vec & sigma,
                      const arma::mat B){
  int N = Z.n_elem;
  arma::mat I_minus_B = - B;
  I_minus_B.diag() += 1; // I - B
  arma::mat Sigma(N,N);
  Simga.diag() = sigma; // Var for noise
  
  // calculating normalizing constant
  // The cov matrix of SAR is : S = (I-B)^{-1} \Sigma (I^T-B^T)^{-1}
  // Thus log(det(S)) = log( det(\Sigma) ) - 2 log (det(I-B))
  double logC = - (N/2) * log (2 * pi) - // 2pi
    .5 * (sum(log(sigma))) + // Sigma
  log(det(I_minus_B)); // I-B part
  
  arma::vec try_solve;
  solve(try_solve, I_minus_B, Z);
  if(!try_solve(0)){
    return(Rcpp::R_NegInf); // if I-B not inversible, return -Inf 
  }
  
  double logH = - 0.5 * X.t() * solve(I_minus_B , Sigma * try_solve); // Hamiltonian
  return(logC+logH);
}

arma::mat SAR_to_Graph(const arma::vec & Bvec,const double & sigma){
  int k = 0.5*(1+sqrt(1+4*Bvec.n_elem));
  arma::mat B = getBmat_helper(Bvec , k);
  arma::mat Sigma(size(B),fill::zeros);
  Sigma.diag() = 1/(sigma*sigma);
  arma::mat I_minus_B = -B;
  I_minus_B.diag() += 1;
  arma::mat full_percision = I_minus_B.t() * Sigma * I_minus_B; //cov matrix of this Sar model
  return(full_percision);
}


Rcpp::List Sigma_to_CAR_Cpp(const arma::mat & Sigma){
  
  arma::sp_mat Q = solve(Sigma); // percision matrix
  arma::mat R = -Sigma; // R matrix, see Ver Hoef et al. 2018
  diag(R) = 0;
  arma::vec M = 1/diag(Q);
  arma::sp_mat invD(size(Q));
  diag(invD) = M ;
  return(Rcpp::List::create(
      Rcpp::Named("C") = invD * R, // C matrix per Ver Hoef et al. 2018
      Rcpp::Named("M") = M // diag of M matrix per Ver Hoef et al. 2018
  ));
}


// tested 20200529
// [[Rcpp::export]]
arma::mat Sample_SAR_Cpp(const int Nsample,
                         const arma::vec & mu,
                         const arma::vec & sigma,
                         const arma::mat B){
  int N = mu.n_elem;
  arma::mat I_minus_B = - B;
  I_minus_B.diag() += 1; // I - B
  arma::mat Sigma(N,N,fill::zeros);
  Sigma.diag() = sigma; // Var for noise
  arma::mat try_solve;
  inv(try_solve , I_minus_B);// This will be all false if I_minus_B is singular
  if(!try_solve(0,0)){
    return(try_solve); // if I-B not inversible, return 0s
  }
  
  arma::mat full_Sigma = try_solve * Sigma * try_solve.t(); //cov matrix of this Sar model
  arma::mat Res;
  Res = mvnrnd(mu, full_Sigma, Nsample ); // this will give column vectors
  return(Res);
}


// inverse Gaussian random variable, tested already 20200528
/*
 * Mitchael,J.R., Schucany, W.R. and Haas, R.W. (1976). Generating
 random roots from variates using transformations with multiple roots.
 American Statistician. 30-2. 88-91.
 */	

// tested already 20200528
// [[Rcpp::export]]
double rinvGau(const double & mu, const double & lambda){
  double b = 0.5 * mu / lambda;
  double a = mu * b;
  double c = 4.0 * mu * lambda;
  double d = mu * mu;
  double res;
  
  if (mu<=0 || lambda<=0) {
    return(NA_REAL);
  }
  
  double u = R::runif(0,1);
  double chisqr = R::rchisq(1);
  
  double x = mu + a * chisqr - b * sqrt( c * chisqr + d * chisqr * chisqr); // solve the smaller root
  res = (u < (mu / (mu + x))) ? x : d/x; // accept the small one with prob mu/(mu+x), otherwise the larger one
  return(res);
}


/*
 * TODO: 
 *  test Gibbs sampler for SAR-lasso
 */

// helper functions


// generate the B function from its vector version, we will not save diagonal since they are alwasy 0
// tested 20200528 
// get convinient X_i: design matrix for each sample
// tested 20200528
arma::sp_mat getDesign_i_helper(const arma::rowvec & X_i,//row vector of that sample
                                int k){ // number of nodes
  int p = X_i.n_elem; // number of predictors
  arma::sp_mat Design(k,k * p);
  for(int j = 0 ; j < k ; ++j){
    Design(j,span( j*p , (j+1) * p - 1 )) = X_i ; 
  }
  return(Design);
}

// tested for dimension compatibility 20200528
arma::vec update_beta_helper(const arma::mat & data,
                             const arma::mat & design,
                             const arma::vec & mu,
                             const double & sigma2,
                             const arma::vec & tau2,
                             const arma::vec & Bvec,
                             int k, int p, int n){
  arma::mat BminusI = getBmat_helper(Bvec,k);
  BminusI.diag() -= 1; // get B-I which is useful
  arma::vec mu_beta(k*p,fill::zeros);
  arma::mat Q_beta(k*p,k*p,fill::zeros);// percision matrix up to sigma^2 scaling
  Q_beta.diag() += 1/tau2;
  arma::vec res;
  
  arma::sp_mat X_i;
  
  /*
   * Q_beta = 1/sigma^2 * (\sum_i [(B-I)X_i]^T[(B-I)X_i] + D_{\tau}^{-1} )
   * 
   * mu_beta = (\sum_i [(B-I)X_i]^T[(B-I)X_i] + D_{\tau}^{-1} )^{-1} (\sum_i [(B-I)X_i]^T[(B-I)(Z_i-\mu)]  )
   * 
   */
  
  
  for(int i = 0 ; i < n ; ++i){
    X_i = getDesign_i_helper(design.row(i),k);
    Q_beta +=  X_i.t() * BminusI.t() * BminusI * X_i ;
    mu_beta +=  X_i.t() * BminusI.t() * (BminusI* (data.col(i)-mu));
  }
  
  mu_beta = solve(Q_beta,mu_beta);
  res = mvnrnd(mu_beta, sigma2 * inv(Q_beta));
  return(res);
}


// tested for dimension 20200528
arma::vec update_B_helper(const arma::mat & data,
                          const arma::mat & design,
                          const arma::vec & mu,
                          const arma::vec & beta,
                          const double & sigma2,
                          const arma::vec & eta2,
                          int k, int n){
  arma::mat Omega(k,k,fill::zeros);
  arma::vec Y_i_tilde;
  arma::vec B(k*(k-1));
  for(int i = 0 ; i < n ; ++i){
    Y_i_tilde = getDesign_i_helper(design.row(i),k) * beta - data.col(i) + mu;
    Omega += Y_i_tilde * Y_i_tilde.t(); // The 
  }
  
  arma::mat Omega_temp;
  arma::mat Sigma_temp;
  arma::mat mu_i;
  arma::uvec ind = linspace<uvec>(0, k-1 , k); // this is a index vector, for deleting the digonal related row and colunms in var cov mat
  arma::uvec ind_remain;
  for(int i = 0 ; i < k ; ++i){
    ind_remain = find(ind != i);
    Omega_temp = Omega.submat(ind_remain,ind_remain); // delete row and column related with the diagnol entry of B
    mu_i = Omega.col(i); // used to solve for expectation
    mu_i = mu_i(ind_remain);
    Omega_temp.diag() += 1/(eta2(span( i*(k-1) , (i+1)*(k-1) - 1 ))); // add LASSO shrik
    
    Sigma_temp = inv(Omega_temp); // calculate cov mat up to sigma^2 scale
    
    mu_i = Sigma_temp * mu_i; // sove for expectation
    
    Sigma_temp = sigma2 * Sigma_temp; // real cov mat
    
    
    B( span(i*(k-1) , (i+1)*(k-1) - 1) ) = mvnrnd(mu_i, Sigma_temp); // save in vector
  }
  
  return(B);
}

// update sigma (standard deviation, not squared)
// test for dimension 20200529

double update_sigma_helper(const arma::mat & data,
                           const arma::mat & design,
                           const arma::vec & mu,
                           const arma::vec & beta,
                           const arma::vec & B,
                           const arma::vec & eta2,
                           const arma::vec & tau2,
                           int k, int p , int n){
  arma::mat Bmat = getBmat_helper(B,k);
  arma::mat betamat = reshape(mat(beta),p,k); // mat version of beta, easier for us to calculate SSE
  arma::mat Xbeta = design * betamat;
  Xbeta = Xbeta.t(); // macro effect of nodes, make sure col as sample
  arma::mat BZ = Bmat * data;
  arma::mat Center_data = data - Xbeta;
  Center_data.each_col() -= mu; // center the macro effect
  arma::mat SSEmat = Center_data - Bmat * Center_data;
  
  double SSE = sum(sum(SSEmat % SSEmat));
  
  double sigma = 1 / sqrt( R::rgamma((n*k-1+k*p+k*(k-1))/2,
                                     2/(SSE+sum(beta % beta / tau2 )+sum(B%B/eta2))) );
  
  return(sigma);
}


//tested dimension 20200529
arma::vec tau2_eta2_update_helper(const arma::vec & beta,
                                  const double & lambda2,
                                  const double & sigma2){
  int n = beta.n_elem;
  arma::vec invtau2(n);
  arma::vec mu_prime = sqrt(lambda2*sigma2/(beta%beta));
  for(int i = 0 ; i < n ; ++i){
    invtau2(i) =  rinvGau(mu_prime(i),lambda2);
  }
  return(1/invtau2);
}

/* 
 * Convention:
 * mcmc obects were vectorized column first, e.g. each row: looks like
 * beta_11,beta_12,beta_13,...,beta_kp
 * 
 */

List SAR_LASSO_Cpp(const arma::mat & data, // raw composition data, column as a sample
                   const arma::mat & design, // design matrix, each ROW as a sample
                   const int n_iter, // how many iteractions?
                   const int n_burn_in, // burn in
                   const int thin_by, // thinning?
                   const double r_beta, // prior on lambda of beta
                   const double delta_beta,
                   const double r_B,
                   const double delta_B,
                   bool progress
){
  int k = data.n_rows; // number of nodes
  int p = design.n_cols; //number of predictors
  int n = data.n_cols; // number of samples
  
  int n_save = floor(n_iter/thin_by); //
  int n_saved = 0;  
  
  // mcmc matrices:
  arma::mat beta_mcmc(n_save,k * p); // beta mcmc
  beta_mcmc += NA_REAL; // initial with NAs
  
  arma::mat B_mcmc(n_save , k * (k - 1)); // column first, but had no diagnol
  B_mcmc += NA_REAL;
  
  arma::mat sigma_mcmc(n_save , 1); //noise standard deviation
  sigma += NA_REAL;
  
  arma::mat mu_mcmc(n_save , k); // mean for node 1 to k
  mu += NA_REAL;
  
  arma::mat lambda(n_save , 2); // LASSO parameter for beta and B
  lambda += NA_REAL;
  
  arma::vec tau2_curr = randg<arma::vec> (k*p,distr_param(r_beta,1/delta_beta)); // current tau^2, for prior of beta
  arma::vec eta2_curr = randg<arma::vec> (k*(k-1),distr_param(r_B,1/delta_B)); // current eta^2, for prior of B
  
  arma::vec mu_curr = mean(data); // current value of mean
  arma::vec beta_curr(k*p , fill::randu); // current value of beta
  arma::vec B_curr(k*(k-1) , fill::randu); // current value of B
  double sigma_curr = sqrt(1/ R::rgamma(1,1)); // current value of sigma
  arma::vec mean_uncertain(k); // for sampling mu
  
  double lambda2_beta = R::rgamma(r_beta,1/delta_beta); // current value of squared LASSO parameter of \beta
  double lambda2_B = R::rgamma(r_B,1/delta_B); // current value of squared LASSO parameter of B
  
  Progress p((n_iter+n_burn_in), progress); // progress bar
  
  for(int i = 0 ; i<(n_iter+n_burn_in) ; ++i){
    if (Progress::check_abort()){
      Rcerr << "keyboard abort\n";
      return(Rcpp::List::create(
          Rcpp::Named("beta") = NA_REAL,
          Rcpp::Named("mu") = NA_REAL,
          Rcpp::Named("B") = NA_REAL,
          Rcpp::Named("sigma") = NA_REAL,
          Rcpp::Named("lambda") = NA_REAL
      );
    }
    // block update start:
    
    // Update betas:
    beta_curr = update_beta_helper(data,design,mu_curr,
                                   sigma_curr * sigma_curr,
                                   tau2_curr,B_curr,
                                   k,p,n);
    
    // Update Bs
    B_curr = update_update_B_helper(data,design,mu_curr,beta_curr,
                                    sigma_curr * sigma_curr,
                                    eta2_curr,
                                    k,n)
      
      // Update mu
      
    mu_curr = mean(data,1) + sqrt(sigma2_curr)/(n*k) * mean_uncertain.randu();
    
    // Update sigma
    // stoped here 20200528
    sigma_curr = update_sigma_helper(data, design, mu_curr,
                                     beta_curr, B_curr,
                                     eta2_curr, tau2_curr,
                                     k,p,n)
      
      // Update tau
    tau2_curr = tau2_eta2_update_helper(beta_curr,lambda2_beta,
                                          sigma_curr * sigma_curr);
    
    
    
    // Update eta
    eta2_curr = tau2_eta2_update_helper(B_curr,lambda2_B,
                                        sigma_curr * sigma_curr);
    
    
    
    // Update lambda_beta
    
    lambda2_beta = R::rgamma(r_beta+k*p,delta_beta+sum(tau2_curr)/2);
    
    
    // Update lambda_B
    lambda2_B = R::rgamma(r_B+k*(k-1),delta_B+sum(eta2_curr)/2);
    
    // saving the state
    if( (i-burn_in)>=0 && (i+1-burn_in)%thin_by ==0 ){
      beta_mcmc.row(n_saved) = beta_curr.t();
      B_mcmc.row(n_saved) = B_curr.t();
      mu_mcmc.row(n_saved) = mu_curr.t();
      sigma.row(n_saved) = sigma_curr.t();
      
      lambda_mcmc(n_save,0) = lambda_beta_curr;
      lambda_mcmc(n_save,1) = lambda_B_curr;
      
      n_save++;
    }
    
    
    p.increment();
  }
  return(Rcpp::List::create(
      Rcpp::Named("beta") = beta_mcmc,
      Rcpp::Named("mu") = mu_mcmc,
      Rcpp::Named("B") = B_mcmc,
      Rcpp::Named("sigma") = sigma_mcmc,
      Rcpp::Named("lambda") = lambda_mcmc
  );
}

