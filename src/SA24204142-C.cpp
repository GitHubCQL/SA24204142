#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp for HW9
//' @param N the number of samples
//' @param a parameter of target density
//' @param b parameter of target density
//' @param n parameter of target density
//' @return a random sample of size \code{N}
//' @examples
//' \dontrun{
//' rnC <- GibbsC(1000, 1, 1, 10)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix GibbsC(int N, double a, double b, int n) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0.2;
  mat(0, 0) = x;
  mat(0, 1) = y;
  for(int i = 1; i < N; i++) {
    x = R::rbinom(n, y);
    y = R::rbeta(x + a, n - x + b);
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}




//' @title Generate samples by Metropolis Adjusted Langevin Algorithm
//' @description This function runs the Metropolis Adjusted Langevin Algorithm provided the \code{logTarget} and gradient \code{glogTarget}.
//' @param N number of samples
//' @param q_init vector of initial values
//' @param sigma standard deviation of proposal distribution
//' @param L number of between-sample random numbers
//' @param logTarget function to calculate the log density of target distribution
//' @param glogTarget function to calculate the gradient of \code{logTarget}
//' @return a list concluding \code{q}, samples of target distribution of size \code{N} and \code{k}, the rejection count 
//' @examples
//' \dontrun{
//' Sigma <- matrix(c(1,0.98,0.98,1), 2, 2)
//' InvSigma <- solve(Sigma)
//' logTarget <- function(x) {
//'   - t(x) %*% InvSigma %*% x / 2
//' }
//' glogTarget <- function(x) {
//'   - InvSigma %*% x
//' }
//' result <- MALA(10000, c(1, 0), 0.18, 20, logTarget, glogTarget)
//' plot(result$q, xlab="x", ylab="y")
//' }
//' @export
// [[Rcpp::export]]
List MALA(int N, NumericVector q_init, double sigma, int L, Function logTarget, Function glogTarget) {
  int d = q_init.size();
  NumericMatrix q(N, d);
  for (int i=0; i<d; i++) {
    q(0, i) = q_init[i];
  }
  int k = 0;
   
  for(int t = 1; t < N; t++){
    NumericVector q0(d);
    for (int i=0; i < d; i++) {
      q0[i] = q(t-1, i);
    }
    NumericVector u = runif(L);
    for (int j=0; j < L; j++) {
      NumericVector mu = q0 + pow(sigma, 2)/2 * NumericVector(glogTarget(q0));
      NumericVector qt = rnorm(d, 0, sigma) + mu;
      NumericVector mut = qt + pow(sigma, 2)/2 * NumericVector(glogTarget(qt));
       
      double proba = exp(as<double>(logTarget(qt))-as<double>(logTarget(q0))
                           -0.5*sum(pow(q0-mut, 2))/pow(sigma, 2)
                           +0.5*sum(pow(qt-mu, 2)/pow(sigma, 2)));
                            
      if(u[j] <= proba){
        q0 = qt;
      }
      else {
        k++;
      }
    }
    for (int i=0; i < d; i++) {
      q(t, i) = q0[i];
    }
  }
   
  List result = List::create(
    Named("q") = q,
    Named("k") = k
  );
  return(result);
}













