// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <climits>
using namespace Rcpp;
using namespace arma;
// FUNCTIONS FOR EXACT SAMPLING //

// Inner function to resize list:
List resize( const List& x, int n ){
    int oldsize = x.size() ;
    List y(n) ;
    for( int i=0; i<oldsize; i++) y[i] = x[i] ;
    return y ;
}

// Inner function to simulate random uniforms in a matrix:
arma::mat RandMat(int nrow, int ncol)
 {
  arma::mat Res = arma::randu(nrow,ncol);
  return(Res);
 }

// Computes maximal and minimal probability of node flipping:
NumericVector PplusMinMax(int i, const arma::sp_mat& J, IntegerVector s, NumericVector h, double beta, IntegerVector responses)
{
  // The function computes the probability that node i is in Response 1 instead of 0, given all other nodes, which might be missing.
  // Output: minimal and maximal probablity
  
  NumericVector H0(2, h[i] * responses[0]); // relevant part of the Hamiltonian for state = 0
  NumericVector H1(2, h[i] * responses[1]); // relevant part of the Hamiltonian for state = 1
  
  NumericVector Res(2);
  
  //int N = J.n_rows;
  NumericVector TwoOpts(2);
  
  for (arma::sp_mat::const_col_iterator it = J.begin_col(i); it != J.end_col(i); ++it)
  {
    if (i != it.row())
    {
      if (s[it.row()] != INT_MIN)
      {
       H0[0] += *it * responses[0] * s[it.row()];
       H0[1] += *it * responses[0] * s[it.row()];
       H1[0] += *it * responses[1] * s[it.row()];
       H1[1] += *it * responses[1] * s[it.row()]; 
      } else 
      {
               
        TwoOpts[0] = *it * responses[1] * responses[0];
        TwoOpts[1] = *it * responses[1] * responses[1];

        if (TwoOpts[1] > TwoOpts[0])
        {
          H1[0] += TwoOpts[0];
          H1[1] += TwoOpts[1];
          
          H0[0] += *it * responses[0] * responses[0];
          H0[1] += *it * responses[0] * responses[1];
        } else 
        {
          H1[0] += TwoOpts[1];
          H1[1] += TwoOpts[0];          
          
          H0[0] += *it * responses[0] * responses[1];
          H0[1] += *it * responses[0] * responses[0];
        }
      }
    }
  }

  Res[0] = exp(beta * H1[0]) / ( exp(beta * H0[0]) + exp(beta * H1[0]) );
  Res[1] = exp(beta * H1[1]) / ( exp(beta * H0[1]) + exp(beta * H1[1]) );
  
  
  return(Res);
}
       
// Inner function:
IntegerVector IsingEx(const arma::sp_mat& graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses, bool exact,
IntegerVector constrain)
{
  // Parameters and results vector:
  int N = graph.n_rows;
  IntegerVector state(N, INT_MIN);
  double u;
  NumericVector P(2);
  int maxChain = 100;
  List U(1);
  int minT = 0;
  bool anyNA = true;
    
  do
  { 
    // Resize U if needed:
    if (minT > 0)
    {
      U = resize(U, minT+1);
    }
    
    // Generate new random numbers:
    U[minT] = RandMat(nIter, N);
    
    // Initialize states:
    for (int i=0; i<N; i++)
    {
      if (exact)
      {
        state[i] = INT_MIN;
      } else 
      {
        state[i] = ifelse(runif(1) < 0.5, responses[1], responses[0])[0];
      }
    }    

    // START ALGORITHM
    for (int t=minT; t > -1;  t--)
    {
      for (int it=0;it<nIter;it++)
      {
        arma::mat Ucur = U[t];
        for (int node=0;node<N;node++)
        {
          u = Ucur(it, node);
          P = PplusMinMax(node, graph, state, thresholds, beta, responses);
          if (u < P[0])
          {
            state[node] = responses[1];
          } else if (u >= P[1])
          {
            state[node] = responses[0];
          } else 
          {
            state[node] = INT_MIN;
          }
        }
      }
    }
    
    anyNA = false;
    if (exact)
    {
      if (minT < maxChain)
      {
       for (int i=0; i<N; i++)
       {
        if (state[i] == INT_MIN)
        {
          anyNA = true;
        }
        } 
      } 
    }    
    minT++;
    
  } while (anyNA);

  // Rf_PrintValue(wrap(minT));
  return(state);
}



// FUNCTIONS FOR METROPOLIS SAMPLER //
double Pplus(int i, const arma::sp_mat& J, IntegerVector s, NumericVector h, double beta, IntegerVector responses)
{
  // The function computes the probability that node i is in Response 1 instead of 0, given all other nodes, which might be missing.
  // Output: minimal and maximal probablity
  
  double H0 = h[i] * responses[0]; // relevant part of the Hamiltonian for state = 0
  double H1 = h[i] * responses[1]; // relevant part of the Hamiltonian for state = 1

  //double Res;

  
  //int N = J.n_rows;
  
  
  for (arma::sp_mat::const_col_iterator it = J.begin_col(i); it != J.end_col(i); ++it)
  {
    if (i != it.row())
    {
       H0 += *it * responses[0] * s[it.row()];
       H1 += *it * responses[1] * s[it.row()];
    }
  }
  
  return(exp(beta * H1) / ( exp(beta * H0) + exp(beta * H1) ));// MH ratio here, need changing
}


IntegerVector IsingMet(const arma::sp_mat& graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses,
IntegerVector constrain)
{
  // Parameters and results vector:
  int N = graph.n_rows;
  IntegerVector state =  ifelse(runif(N) < 0.5, responses[1], responses[0]);
  for (int i=0; i<N; i++)
  {
    if (constrain[i] != INT_MIN)
    {
      state[i] = constrain[i];
    }
  }
  double u;
  double P;
    
    // START ALGORITHM
    for (int it=0;it<nIter;it++)
    {
      for (int node=0;node<N;node++)
      {
        if (constrain[node] == INT_MIN)
        {
         u = runif(1)[0];
         P = Pplus(node, graph, state, thresholds, beta, responses);
          if (u < P)
         {
           state[node] = responses[1];
         } else 
         {
           state[node] = responses[0];
         } 
        }
      }
    }
   
  return(state);
}


///ISING PROCESS SAMPLER:
// [[Rcpp::export]]
IntegerMatrix IsingProcess(int nSample, const arma::sp_mat& graph, NumericVector thresholds, double beta, IntegerVector responses)
{
  // Parameters and results vector:
  int N = graph.n_rows;
  IntegerVector state =  ifelse(runif(N) < 0.5, responses[1], responses[0]);
  double u;
  double P;
  IntegerMatrix Res(nSample,N);
  int node;
    
    // START ALGORITHM
    for (int it=0;it<nSample;it++)
    {
      node = floor(R::runif(0,N));
        u = runif(1)[0];
        P = Pplus(node, graph, state, thresholds, beta, responses);
        if (u < P)
        {
          state[node] = responses[1];
        } else 
        {
          state[node] = responses[0];
        }
        for (int k=0; k<N; k++) Res(it,k) = state[k];
    }
   
  return(Res);
}





// OVERAL FUNCTION //
// [[Rcpp::export]]
IntegerMatrix IsingSamplerCpp(int n, const arma::sp_mat& graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses, bool exact,
IntegerVector constrain)
{
  int Ni = graph.n_rows;
  IntegerMatrix Res(n,Ni);
  IntegerVector state(Ni);
  //IntegerVector constrainVec(Ni);
  if (exact)
  {
    for (int s=0;s<n;s++)
    {
      //for (int i=0;i<Ni;i++) constrainVec[i] = constrain(s,i);
      state = IsingEx(graph, thresholds, beta, nIter, responses, exact, constrain);
      for (int i=0;i<Ni;i++) Res(s,i) = state[i];
    }
  } else 
  {
    for (int s=0;s<n;s++)
    {
      //for (int i=0;i<Ni;i++) constrainVec[i] = constrain(s,i);
      state = IsingMet(graph, thresholds, beta, nIter, responses, constrain);
      for (int i=0;i<Ni;i++) Res(s,i) = state[i];
    }
  }
  
  return(Res);
}


// HELPER FUNCTIONS //
// Hamiltonian:
// [[Rcpp::export]]
double H(const arma::sp_mat& J, IntegerVector s, NumericVector h)
{
  double Res = 0;
  int N = J.n_rows;
  for (int i=0;i<N;i++)
  {
    Res -= h[i] * s[i];
    for (arma::sp_mat::const_col_iterator it = J.begin_col(i); it != J.end_col(i); ++it)
    {
      if (it.row()!=i) Res -= *it * s[i] * s[it.row()] * .5;
    }
  }
  return(Res);
}


// [[Rcpp::export]]
double PartitionCpp(
    const arma::sp_mat & graph,
    const NumericVector & thr,
    const double & beta,
    const IntegerVector &responses){
  
  int N = graph.n_rows;
  int n_possible = pow(2,N);
  int t=0;
  double Z=0;
  
  IntegerVector temp(N);
  
  for(int i = 0; i<n_possible ; ++i){
    t = i;
    for(int j = 0; j<N ; ++j){
      temp[j] = responses[t % 2];//use binary number coding
      t = t>>1;
    }
    Z+=exp(-beta * H(graph,temp,thr));
  }
  
  return(Z);
  
  
}

// calculating state probability directly 
// [[Rcpp::export]]
double IsingStateProbCpp(const arma::vec & Z, 
                         const arma::sp_mat & graph,
                         const arma::vec & thr,
                         const IntegerVector &responses){
  double partion = PartitionCpp(graph,
                          Rcpp::NumericVector(thr.begin(),
                                              thr.end()),
                          1,responses);
  return(exp(- H(graph, IntegerVector(Z.begin(),
                                     Z.end()),                       
                       Rcpp::NumericVector(thr.begin(),
                                              thr.end())))/partion);
} // passed 2/4/2020



// [[Rcpp::export]]
double ColSumsarma(arma::mat A){
  return(sum( is_na(as<NumericMatrix>(wrap(A)))));
}


// Likelihood without Z
// [[Rcpp::export]]
double f(IntegerMatrix Y, const arma::sp_mat& J, NumericVector h)
{
  double Res = 1;
  int Np = Y.nrow();
  int Ni = J.n_cols;
  IntegerVector s(Ni);
  for (int p=0;p<Np;p++)
  {
    for (int i=0;i<Ni;i++)
    {
      s[i] = Y(p,i);
    }
    Res *= exp(-1.0 * H(J, s, h));
  }
  return(Res);
}

// [[Rcpp::export]]
double Ising_PseudoLikelihood_Cpp(const arma::mat & x, 
                                  const arma::sp_mat & graph, 
                                  const arma::vec & thresholds, 
                                   
                                  const IntegerVector & responses, 
                                  bool log){
  int n = x.n_rows;
  int k = x.n_cols;
  
  double logPS = 0;
  double e = 0;
  
  // for every response, compute log PL and sum:
  for (int i = 0 ; i < n ; ++i){
    // for every node, compute log PL and sum:
    for (int j = 0 ; j < k ; ++j){
      e = thresholds[j] + 0.5 * x.row(i) * graph.col(j);
      logPS += x(i,j) * e  -  log(exp(responses[0]*e) + exp(responses[1]*e));
    }
  }
  
  
  return(logis? logPS : exp(logPS));
}



