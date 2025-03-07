#' Estimate Parameters for a Single Run of the Factor OUSSM Model
#'
#' This function estimates the parameters of the Factor Ornstein-Uhlenbeck 
#' State-Space Model (OUSSM) for a single dataset using Maximum Likelihood Estimation (MLE).
#'
#' @param theta.off.diag.init A numeric vector representing the initial values 
#'   for the off-diagonal elements of the state transition matrix \code{Theta}.
#' @param theta.diag.init A numeric vector representing the initial values for 
#'   the diagonal elements of \code{Theta}.
#' @param ytn.sim A matrix of dimensions \code{(N, p.dim)} containing 
#'   the observed data for estimation.
#' @param N An integer specifying the number of time points.
#' @param p.dim An integer specifying the number of observed variables.
#' @param m.dim An integer specifying the number of latent state variables.
#' @param tn.diff A numeric vector of length \code{N-1} representing the time 
#'   differences between consecutive observations.
#' @return A list containing:
#' \item{theta.matrix.dist}{The estimated \code{Theta} matrix.}
#' \item{C.matrix.dist}{The estimated state transition matrix \code{C}.}
#' \item{Z.matrix.dist}{The estimated observation matrix \code{Z}.}
#' \item{H.matrix.dist}{The estimated observation noise covariance matrix \code{H}.}
#' \item{mu.dist}{The estimated mean vector.}
#' \item{Hessian.record}{The Hessian matrix from optimization.}
#' \item{L.dist}{The log-likelihood value from optimization.}
#' \item{result.para.dist}{The estimated parameter vector.}
#' \item{AIC}{The Akaike Information Criterion (AIC) value for model selection.}
#' \item{BIC}{The Bayesian Information Criterion (BIC) value for model selection.}
#' @details
#' The function follows these steps:
#' \itemize{
#'   \item Computes initial estimates for \code{Theta}, \code{Z}, and \code{H} based on the covariance of the data.
#'   \item Uses eigen decomposition to initialize the observation matrix \code{Z}.
#'   \item Runs the \code{optim()} function with the \code{L-BFGS-B} method to estimate parameters.
#'   \item Converts the estimated parameters into structured matrices using \code{para.to.model()}.
#'   \item Ensures sign consistency using \code{convert.Z.theta()}.
#'   \item Computes model selection criteria (AIC and BIC).
#' }
#' @examples
#' # Example usage (assuming required parameters and data are predefined)
#' result <- estimation.OU.single(theta.off.diag.init, theta.diag.init, ytn.sim, 
#'                                N, p.dim, m.dim, tn.diff)
#' result$theta.matrix.dist  # Access the estimated Theta matrix
#' @export
estimation.OU.single <- function(theta.off.diag.init, theta.diag.init, 
                                 ytn.sim, N, p.dim, m.dim, tn.diff) {
  # Convert ytn.sim to a matrix format
  temp.sim <- as.matrix(ytn.sim)
  
  # Compute covariance-based initialization
  M <- cov(temp.sim[1:(N-1), ], temp.sim[2:N, ])
  M <- eigen(M + t(M))
  
  # Initialize Z matrix based on eigenvectors
  if (p.dim >= m.dim) {
    Z.init <- M$vectors[, 1:m.dim]
  } else {
    Z.init <- cbind(M$vectors, matrix(1, nrow = p.dim, ncol = (m.dim - p.dim)))
  }
  Z.init <- c(t(Z.init))
  
  # Initialize other parameters
  H.par.init <- diag(cov(temp.sim)) / 3
  mu.init <- colMeans(temp.sim)
  
  # Combine into initial parameter vector
  init.par <- c(theta.off.diag.init, log(theta.diag.init), Z.init, H.par.init, mu.init)
  
  # Run MLE optimization using L-BFGS-B
  result.sim <- optim(init.par, p.dim = p.dim, m.dim = m.dim, N = N, 
                      ytn = temp.sim, tn.diff = tn.diff, neglogL, 
                      method = "L-BFGS-B", hessian = TRUE)
  
  # Convert estimated parameters into structured matrices
  result.mat <- para.to.model(result.sim$par, m.dim, p.dim)
  
  # Ensure sign consistency
  convert.result <- convert.Z.theta(result.mat$theta, result.mat$Z.mat)
  
  # Store estimated matrices
  theta.matrix.dist <- convert.result$theta
  C.matrix.dist <- convert.result$C.matrix
  Z.matrix.dist <- convert.result$Z.matrix
  H.matrix.dist <- result.mat$H.mat
  mu.dist <- result.mat$mu
  Hessian.record <- result.sim$hessian
  L.dist <- result.sim$value
  result.para.dist <- result.sim$par
  
  # Compute model selection criteria
  AIC.value <- 2 * L.dist + 2 * length(result.para.dist)
  BIC.value <- 2 * L.dist + log(N) * length(result.para.dist)
  
  # Return results
  return(list("theta.matrix.dist" = theta.matrix.dist, 
              "C.matrix.dist" = C.matrix.dist, 
              "Z.matrix.dist" = Z.matrix.dist, 
              "H.matrix.dist" = H.matrix.dist, 
              "mu.dist" = mu.dist, 
              "Hessian.record" = Hessian.record,  
              "L.dist" = L.dist, 
              "result.para.dist" = result.para.dist, 
              "AIC" = AIC.value, 
              "BIC" = BIC.value))
}
