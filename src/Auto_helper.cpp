// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <tgmath.h>
using namespace Rcpp;
static const double pi = 3.141592653589793238462643383280;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

/* Convention:
 *    single sample will be saved as col vectors
 */

// Auto-Poisson

arma::vec Auto_Poisson_Gibbs_Single_Cpp(const arma::sp_mat &graph, // graph
                                        const arma::vec &thresholds, // mean 
                                        const int & Winsorized, // right censor
                                        const int & nIter){  // number of Gibbs iters
  int N = graph.n_rows;
  int temp;
  arma::vec Res(N);
  
  for (int i = 0 ; i < N ; ++i){
    Res(i) = R::rpois(thresholds(i)); // initializing 
  }
  
  for(int i = 0 ; i < nIter ; ++i){ // outer loop, for number of samples
    for(int j = 0 ; j < N ; ++j){
      temp = R::rpois(exp(thresholds(j) + 0.5 * (Res.t() * graph.col(j)))); // Gibbs scan
      Res(i) = temp > Winsorized ? Winsorized : temp; // right censored 
    }
  }
  return(Res);
} 


// [[Rcpp::export]]
arma::mat Auto_Poisson_Gibbs_Batch_Cpp(const int & Nsample, // number of samples
                                       const arma::sp_mat &graph, // graph
                                       const arma::vec &thresholds, // mean 
                                       const int & Winsorized, // right censor
                                       const int & nIter){  
  int N = graph.n_rows;
  arma::mat Res(N,Nsample)
  for(int i = 0 ; i < Nsample ; ++i ){
    Res.col(i) = Auto_Poisson_Gibbs_Single_Cpp(graph,thresholds,Winsorized,nIter);
  }
  return(Res);
}

// [[Rcpp::export]]
double Auto_Poisson_Pseudo_likelihood_Cpp(const arma::mat Sample,
                                          const arma::sp_mat & graph,
                                          const arma::vec & thresholds){
  int Nsample = Sample.n_cols; //samples
  int N = Sample.n_rows; // nodes
  double Res = 0;
  arma::mat Interactions = graph * Sample ; // interaction terms
  for(int i = 0 ; i < Nsample ; ++i){
    arma::vec Lambda = exp(thresholds + Interactions.col(i)); // lambda of Poisson of that repeat
    for(int j = 0 ; j < N ; ++j){
      Res += R::dpois(Sample(j,i),Lambda(j),true); // likelihood, use Poisson approximation 
    } 
  }
  
  return(Res);
}

// Q-function for auto poisson 
double Auto_Poisson_Q_Cpp(const arma::mat Sample,
                          const arma::sp_mat & graph,
                          const arma::vec & thresholds){
  // not sure if this line work sonce sum(...) will give a row vec
  arma::rowvec thr_contribution = thresholds.t() * Sample - sum( lgamma(Sample + 1),0); // non-interaction part of Besag's Q-function, see Augustin et al. 2006
  double Res = 0;
  int Nsample = Sample.n_cols;
  
  for(int i = 0 ; i < Nsample ; ++i){
    Res += thr_contribution(i) + .5 * (Sample.col(i).t() * graph * Sample.col(i)); // Besag's Q function for (spatial) auto Poisson model 
  }
  return(Res);
}


// logLik of a normal response using composition as predictor
double logLik_Normal_Response(const arma::vec & data,
                              const arma::mat & Design,
                              const arma::vec & beta_fix,
                              const arma::mat & composition,
                              const arma::vec & beta_composition,
                              const double & sigma){
  int N = Design.n_rows;
  arma::vec mu = Design * beta_fix + composition.t() * beta_composition;
  double Res = 0;
  for(int i = 0 ; i < N ; ++i){
    Res += R::dnorm(data(i) , mu(i) , true);
  }
  return(Res);
}


// Laplace distribution density, for L1 regularization
// [[Rcpp::export]]
double dLaplace_vec_Cpp(const arma::vec &x,const double &mu, const double &b, bool log){
  double logd = sum(- abs(x-mu)/b - log(2*b));
  return(log ? logd : exp(logd));
}

// [[Rcpp::export]]
double dLaplace_Cpp(const double &x,const double &mu, const double &a, bool log){
  double logd = - a * abs(x-mu) + log(a/2);
  return(log ? logd : exp(logd));
}


/* TODO:
 * 1. Coupling from the past for one perfect sample from the auto poisson model
 * 2. Auto-NegBin model sampling, 
 * 3. derive the Q-function 
 */



// Auto-normal

/* Starting here is auto-normal, the main idea was fitting a SAR and convert it to CAR if needed, see Ver Hoef et al. 2018
 * The Reason for using SAR rather than CAR is that SAR asks for a ralative easier to check condition 
 *   condition on B (I-B) non-singular, while CAR asks for C to have positive eigen value
 *   We can use independent Laplace prior in SAR to get similar result as LASSO.
 * This is all because Sigma has to be positive semi-definite.
 * Also, matrix B can be asymmertic, 
 *   however this will cause some explaning problem, 
 *    as B is no longer an adjacency matrix of the faithfuil graph related with the joint distribution 
 */


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


// [[Rcpp::export]]
arma::mat Sample_SAR_Cpp(const int Nsample,
                         const arma::vec & mu,
                           const arma::vec & sigma,
                         const arma::mat B){
  int N = mu.n_elem;
  arma::mat I_minus_B = - B;
  I_minus_B.diag() += 1; // I - B
  arma::mat Sigma(N,N);
  Simga.diag() = sigma; // Var for noise
  arma::mat try_solve;
    solve(try_solve, I_minus_B);// This will be all false if I_minus_B is singular
  if(!try_solve(0,0)){
    return(try_solve); // if I-B not inversible, return 0s
  }
  
  arma::mat full_Sigma = try_solve * Sigma * try_solve.t(); //cov matrix of this Sar model
  arma::mat Res;
  mvnrnd(Res, mu, full_Sigma, Nsample )); // this will give column vectors
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
 *  Gibbs sampler for auto guassian lasso
 */

// helper functions


// generate the B function from its vector version, we will not save diagonal since they are alwasy 0
// tested 20200528 
arma::mat gettingBmat_helper(const arma::vec &B, int k){
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
  arma::mat BminusI = gettingBmat_helper(Bvec,k);
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




/* 
 * Convention:
 * mcmc obects were vectorized column first, e.g. each row: looks like
 * beta_11,beta_12,beta_13,...,beta_kp
 * 
 */

List Auto_Gaussian_LASSO_sampler_Cpp(const arma::mat & data, // raw composition data, column as a sample
                                     const arma::mat & design, // design matrix, each ROW as a sample
                                     const int n_iter, // how many iteractions?
                                     const int n_burn_in, // burn in
                                     const int thin_by, // thinning?
                                     const double r_beta, // prior on lambda of beta
                                     const double delta_beta,
                                     const double r_B,
                                     const double delta_B
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
  
  arma::mat sigma(n_save , 1); //noise standard deviation
  sigma += NA_REAL;
  
  arma::mat mu(n_save , k); // mean for node 1 to k
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
  
  Progress p((n_iter+n_burn_in), display_progress); // progress bar
  
  for(int i = 0 ; i<(n_iter+n_burn_in) ; ++i){
    if (Progress::check_abort()) return NA_REAL;
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
    
    mu = mean(data,1) + sqrt(sigma2_curr)/(n*k) * mean_uncertain.randu();
    
    // Update sigma
    
    
    // Update tau
    
    
    // Update eta
    
    
    // Update lambda_beta
    
    
    // Update lambda_B
    
    
    // saving the state
    if( (i-burn_in)>=0 && (i+1-burn_in)%thin_by ==0 ){
      beta_mcmc.row(n_saved) = beta_curr;
      B_mcmc.row(n_saved) = B_curr;
      mu_mcmc.row(n_saved) = mu_curr;
      sigma.row(n_saved) = sigma_curr;
      
      lambda_mcmc(n_save,0) = lambda_beta_curr;
      lambda_mcmc(n_save,1) = lambda_B_curr;
      
      n_save++;
    }
    
    
    p.increment();
  }
  
}





