
#' @importFrom microbenchmark microbenchmark
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_point geom_smooth geom_raster geom_vline stat_function labs theme_bw theme scale_color_manual scale_fill_gradientn element_blank element_text xlim ylim
#' @importFrom RColorBrewer brewer.pal
#' @import patchwork
#' @importFrom MASS mvrnorm
#' @import boot
#' @import bootstrap
#' @import DAAG
#' @importFrom lpSolve lp
#' @import Rcpp
#' @importFrom graphics abline par lines
#' @importFrom stats rnorm dnorm qnorm runif rpois rgamma dgamma rexp dexp rbeta pbeta dcauchy qcauchy rbinom lm kmeans ecdf cov var quantile cor.test t.test chisq.test p.adjust integrate uniroot qqplot
#' @useDynLib SA24204142
NULL


#' @title Generate samples by Random-Walk Metropolis Algorithm 
#' @description This function runs the Random Walk Metropolis Algorithm provided the \code{logTarget}.
#' @param N number of samples
#' @param q_init vector of initial values
#' @param sigma standard deviation of proposal distribution
#' @param L number of between-sample random numbers
#' @param logTarget function to calculate the log density of target distribution
#' @param isRandom logical value to determine whether to select parameter \code{sigma} randomly
#' @param interval interval from which we select \code{sigma} randomly
#' @return a list concluding \code{q}, samples of target distribution of size \code{N} and \code{k}, the rejection count 
#' @examples
#' \dontrun{
#'   Sigma <- matrix(c(1,0.98,0.98,1), 2, 2)
#'   InvSigma <- solve(Sigma)
#'   logTarget <- function(x) {
#'     - t(x) %*% InvSigma %*% x / 2  
#'   }
#'   result <- MH(400, c(1, 0), 0.18, 20, logTarget)
#'   plot(result$q, xlab="x", ylab="y")
#' }
#' @export
MH <- function(N, q_init, sigma=0.1, L=1, logTarget, isRandom=FALSE, interval=c(0, 1)){
  d <- length(q_init)
  q <- matrix(nrow=N, ncol=d)
  q[1, ] <- q_init
  k <- 0
  for(i in 2:N){
    q0 <- q[i-1, ]
    u <- runif(L)
    for (j in 1:L) {
      if(isRandom) {
        sd <- runif(1, interval[1], interval[2])
      }
      else{
        sd <- sigma
      }
      y <- rnorm(d, q0, sd)
      proba <- exp(logTarget(y)-logTarget(q0))
      if(u[j] <= proba){
        q0 <- y
      }
      else{
        k <- k+1
      }
    }
    q[i, ] <- q0
  }
  return(list(q=q, k=k))
}



LeapFrog <- function(q0, p0, delta, L, glogTarget) {
  q <- q0
  p <- p0
  for (i in 1:L) {
    p <- p + 0.5 * delta * glogTarget(q)
    q <- q + delta * p
    p <- p + 0.5 * delta * glogTarget(q) 
  }
  return(list(q=q, p=p))
}



#' @title Generate samples by Hamiltonian Monte Carlo
#' @description This function runs the HMC algorithm provided the \code{logTarget} and gradient \code{glogTarget}.
#' @param N number of samples
#' @param q_init vector of initial values
#' @param delta step size for \code{leapfrog}
#' @param L number of \code{leapfrog} steps in each iteration
#' @param logTarget function to calculate the log density of target distribution
#' @param glogTarget function to calculate the gradient of \code{logTarget}
#' @param isRandom logical value to determine whether to select parameter \code{delta} randomly
#' @param interval interval from which we select \code{delta} randomly
#' @return a list concluding \code{q}, samples of target distribution of size \code{N} and \code{k}, the rejection count 
#' @examples
#' \dontrun{
#'   Sigma <- matrix(c(1,0.98,0.98,1), 2, 2)
#'   InvSigma <- solve(Sigma)
#'   logTarget <- function(x) {
#'     - t(x) %*% InvSigma %*% x / 2
#'   }
#'   glogTarget <- function(x) {
#'     - InvSigma %*% x
#'   }
#'   result <- HMC(10000, c(1, 0), 0.18, 20, logTarget, glogTarget)
#'   plot(result$q, xlab="x", ylab="y")
#' }
#' @export
HMC <- function(N, 
                q_init, 
                delta=0.1, 
                L=10, 
                logTarget, 
                glogTarget, 
                isRandom=FALSE, 
                interval=c(0, 1)){
  k = 0
  d <- length(q_init)
  q <- matrix(nrow=N, ncol=d)
  q[1, ] <- q_init
  
  for(t in 2:N){
    q0 <- q[t-1, ]
    p0 <- rnorm(d)
    if (isRandom) {
      epsilon <- runif(1, interval[1], interval[2])
    }
    else {
      epsilon <- delta
    }
    tmp <- LeapFrog(q0, p0, epsilon, L, glogTarget)
    qt <- tmp$q
    pt <- tmp$p
    proba <- exp(logTarget(qt)-logTarget(q0)-0.5*sum(pt^2)+0.5*sum(p0^2))
    u <- runif(1)
    if(u <= proba) {
      q[t, ] <- qt
    }
    else {
      q[t, ] <- q[t-1, ]
      k <- k+1
    }
  }
  return(list(q=q, k=k))
}



BuildTree <- function(q, p, u, v, j, delta, logTarget, glogTarget) {
  if (j == 0) {
    tmp <- LeapFrog(q, p, v*delta, 1, glogTarget)
    q1 <- tmp$q
    p1 <- tmp$p
    n1 <- as.numeric(log(u)<=(logTarget(q1)-0.5*sum(p1^2)))
    s1 <- as.numeric((logTarget(q1)-0.5*sum(p1^2)) > (log(u)-1000))
    return(list(q_minus=q1, 
                p_minus=p1, 
                q_plus=q1, 
                p_plus=p1, 
                q_prime=q1, 
                n=n1, 
                s=s1))
  }
  
  else {
    tmp <- BuildTree(q, p, u, v, j-1, delta, logTarget, glogTarget)
    q_minus <- tmp$q_minus
    p_minus <- tmp$p_minus
    q_plus <- tmp$q_plus
    p_plus <- tmp$p_plus
    q1 <- tmp$q_prime
    n1 <- tmp$n
    s1 <- tmp$s
    
    if (s1 == 1) {
      if (v == -1) {
        tmp <- BuildTree(q_minus, p_minus, u, v, j-1, delta, logTarget, glogTarget)
        q_minus <- tmp$q_minus
        p_minus <- tmp$p_minus
        q2 <- tmp$q_prime
        n2 <- tmp$n
        s2 <- tmp$s
      }
      
      else {
        tmp <- BuildTree(q_plus, p_plus, u, v, j-1, delta, logTarget, glogTarget)
        q_plus <- tmp$q_plus
        p_plus <- tmp$p_plus
        q2 <- tmp$q_prime
        n2 <- tmp$n
        s2 <- tmp$s
      }
      
      proba <- n2 / max(n1 + n2, 1)
      if (runif(1) <= proba) {
        q1 <- q2
      }
      s1 <- s2 * (sum((q_plus-q_minus)*p_minus) >= 0)*(sum((q_plus-q_minus)*p_plus) >= 0)
      n1 <- n1 + n2
    }
    
    return(list(q_minus=q_minus, 
                p_minus=p_minus, 
                q_plus=q_plus, 
                p_plus=p_plus, 
                q_prime=q1, 
                n=n1, 
                s=s1))
  }
}



#' @title Generate samples by No U-Turn Sampler
#' @description This function runs the NUTS algorithm provided the \code{logTarget} and gradient \code{glogTarget}.
#' @param N number of samples
#' @param q_init vector of initial values
#' @param delta step size for \code{leapfrog}
#' @param logTarget function to calculate the log density of target distribution
#' @param glogTarget function to calculate the gradient of \code{logTarget}
#' @return a matrix concluding samples of target distribution of size \code{N} 
#' @examples
#' \dontrun{
#'   Sigma <- matrix(c(1,0.98,0.98,1), 2, 2)
#'   InvSigma <- solve(Sigma)
#'   logTarget <- function(x) {
#'     - t(x) %*% InvSigma %*% x / 2
#'   }
#'   glogTarget <- function(x) {
#'     - InvSigma %*% x
#'   }
#'   result <- NUTS(10000, c(1, 0), 0.18, logTarget, glogTarget)
#'   plot(result, xlab="x", ylab="y")
#' }
#' @export
NUTS <- function(N, q_init, delta, logTarget, glogTarget) {
  d <- length(q_init)

  qs <- matrix(nrow=N, ncol=d)
  qs[1, ] <- q_init  
  
  for (t in 2:N) {
    p0 <- rnorm(d)
    u <- runif(1, 0, exp(logTarget(qs[t-1,])-0.5*sum(p0^2)))
    q_minus <- qs[t-1, ]
    q_plus <- qs[t-1, ]
    p_minus <- p0
    p_plus <- p0
    j <- 0
    qs[t, ] <- qs[t-1, ]
    n <- 1
    s <- 1
    
    while (s == 1) {
      v <- sample(c(-1, 1), 1)
      if (v == -1) {
        tmp <- BuildTree(q_minus, p_minus, u, v, j, delta, logTarget, glogTarget)
        q_minus <- tmp$q_minus
        p_minus <- tmp$p_minus
        q1 <- tmp$q_prime
        n1 <- tmp$n
        s1 <- tmp$s
      }
      
      else {
        tmp <- BuildTree(q_plus, p_plus, u, v, j, delta, logTarget, glogTarget)
        q_plus <- tmp$q_plus
        p_plus <- tmp$p_plus
        q1 <- tmp$q_prime
        n1 <- tmp$n
        s1 <- tmp$s
      }
      
      if (s1 == 1) {
        if (runif(1) <= n1/n) {
          qs[t, ] <- q1
        }
      }
      n <- n + n1
      s <- s1 * (sum((q_plus-q_minus)*p_minus) >= 0) * (sum((q_plus-q_minus)*p_plus) >= 0)
      j <- j + 1
    }
  }
  return(qs)
}

