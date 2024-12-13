# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A Gibbs sampler using Rcpp
#' @description A Gibbs sampler using Rcpp for HW9
#' @param N the number of samples
#' @param a parameter of target density
#' @param b parameter of target density
#' @param n parameter of target density
#' @return a random sample of size \code{N}
#' @examples
#' \dontrun{
#' rnC <- GibbsC(1000, 1, 1, 10)
#' par(mfrow=c(2,1));
#' plot(rnC[,1],type='l')
#' plot(rnC[,2],type='l')
#' }
#' @export
GibbsC <- function(N, a, b, n) {
    .Call('_SA24204142_GibbsC', PACKAGE = 'SA24204142', N, a, b, n)
}

#' @title Generate samples by Metropolis Adjusted Langevin Algorithm
#' @description This function runs the Metropolis Adjusted Langevin Algorithm provided the \code{logTarget} and gradient \code{glogTarget}.
#' @param N number of samples
#' @param q_init vector of initial values
#' @param sigma standard deviation of proposal distribution
#' @param L number of between-sample random numbers
#' @param logTarget function to calculate the log density of target distribution
#' @param glogTarget function to calculate the gradient of \code{logTarget}
#' @return a list concluding \code{q}, samples of target distribution of size \code{N} and \code{k}, the rejection count 
#' @examples
#' \dontrun{
#' Sigma <- matrix(c(1,0.98,0.98,1), 2, 2)
#' InvSigma <- solve(Sigma)
#' logTarget <- function(x) {
#'   - t(x) %*% InvSigma %*% x / 2
#' }
#' glogTarget <- function(x) {
#'   - InvSigma %*% x
#' }
#' result <- MALA(10000, c(1, 0), 0.18, 20, logTarget, glogTarget)
#' plot(result$q, xlab="x", ylab="y")
#' }
#' @export
MALA <- function(N, q_init, sigma, L, logTarget, glogTarget) {
    .Call('_SA24204142_MALA', PACKAGE = 'SA24204142', N, q_init, sigma, L, logTarget, glogTarget)
}

