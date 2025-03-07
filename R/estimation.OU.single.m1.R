#' Estimate Parameters for One-Dimensional OUSSM Using Maximum Likelihood
#'
#' This function estimates the parameters of a one-dimensional Ornstein-Uhlenbeck 
#' state-space model (OUSSM) using maximum likelihood estimation (MLE) via 
#' the Kalman filter and the L-BFGS-B optimization method.
#'
#' @param theta.init Initial value for the mean-reverting rate (\eqn{\theta}).
#' @param ytn.sim A numeric matrix of observed time-series data (size: \code{N} x \code{p.dim}).
#' @param N An integer representing the number of time points in the dataset.
#' @param p.dim An integer specifying the number of observed variables.
#' @param tn.diff A numeric vector of time differences between observations.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{theta.matrix.dist}: Estimated mean-reverting rate (\eqn{\theta}).
#'   \item \code{C.matrix.dist}: Estimated state transition matrix.
#'   \item \code{Z.matrix.dist}: Estimated loading matrix (adjusted for sign consistency).
#'   \item \code{H.matrix.dist}: Estimated measurement noise covariance matrix.
#'   \item \code{mu.dist}: Estimated mean vector.
#'   \item \code{Hessian.record}: Hessian matrix from optimization.
#'   \item \code{L.dist}: Negative log-likelihood value.
#'   \item \code{result.para.dist}: Vector of estimated parameters.
#'   \item \code{AIC}: Akaike Information Criterion (AIC) value.
#'   \item \code{BIC}: Bayesian Information Criterion (BIC) value.
#' }
#'
#' @details
#' The function follows these steps:
#' \itemize{
#'   \item Computes initial estimates of \code{Z.init}, \code{H.par.init}, and \code{mu.init}.
#'   \item Optimizes the negative log-likelihood function using \code{optim()} with L-BFGS-B.
#'   \item Converts the estimated parameter vector into structured matrices.
#'   \item Adjusts the sign of the estimated loading matrix (\code{Z.matrix.dist}) for consistency.
#'   \item Computes AIC and BIC values for model selection.
#' }
#'
#' @examples
#' # Example usage (assuming required inputs are predefined)
#' theta.init <- 0.5
#' ytn.sim <- matrix(rnorm(100), nrow = 50, ncol = 2)
#' N <- 50
#' p.dim <- 2
#' tn.diff <- rep(1, N - 1)
#' result <- estimation.OU.single.m1(theta.init, ytn.sim, N, p.dim, tn.diff)
#'
#' @export
estimation.OU.single.m1 <- function(theta.init, ytn.sim, N, p.dim, tn.diff) {
  
  # Convert data to matrix format
  temp.sim <- as.matrix(ytn.sim)
  
  # Compute initial estimates
  M <- cov(temp.sim[1:(N - 1),], temp.sim[2:N,])
  M <- eigen(M + t(M))
  Z.init <- M$vectors[, 1]
  Z.init <- c(t(Z.init))
  H.par.init <- diag(cov(temp.sim)) / 3
  mu.init <- colMeans(temp.sim)
  
  # Construct initial parameter vector
  init.par <- c(log(theta.init), Z.init, H.par.init, mu.init)
  
  # Perform maximum likelihood estimation using L-BFGS-B optimization
  result.sim <- optim(init.par, 
                      p.dim = p.dim, 
                      N = N, 
                      ytn = temp.sim, 
                      tn.diff = tn.diff, 
                      neglogL.m1, 
                      method = "L-BFGS-B", 
                      hessian = TRUE)
  
  # Convert estimated parameters into model matrices
  result.mat <- para.to.model.m1(result.sim$par, p.dim)
  
  # Adjust sign consistency in the loading matrix
  theta.matrix.dist <- result.mat$theta
  C.matrix.dist <- result.mat$C.matrix
  Z.matrix.dist <- result.mat$Z.mat[1,] / abs(result.mat$Z.mat[1,]) * result.mat$Z.mat
  H.matrix.dist <- result.mat$H.mat
  mu.dist <- result.mat$mu
  
  # Store optimization results
  Hessian.record <- result.sim$hessian
  L.dist <- result.sim$value
  result.para.dist <- result.sim$par
  
  # Compute AIC and BIC values
  AIC.value <- 2 * L.dist + 2 * length(result.para.dist)
  BIC.value <- 2 * L.dist + log(N) * length(result.para.dist)
  
  # Return structured output
  result <- list("theta.matrix.dist" = theta.matrix.dist, 
                 "C.matrix.dist" = C.matrix.dist, 
                 "Z.matrix.dist" = Z.matrix.dist, 
                 "H.matrix.dist" = H.matrix.dist, 
                 "mu.dist" = mu.dist, 
                 "Hessian.record" = Hessian.record, 
                 "L.dist" = L.dist, 
                 "result.para.dist" = result.para.dist, 
                 "AIC" = AIC.value, 
                 "BIC" = BIC.value)
  
  return(result)
}
