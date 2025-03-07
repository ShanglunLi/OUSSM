#' Estimate Parameters for the Factor OUSSM Model
#'
#' This function estimates the parameters of the Factor Ornstein-Uhlenbeck 
#' State-Space Model (OUSSM) using Maximum Likelihood Estimation (MLE). It 
#' initializes parameters, runs the optimization using the Kalman filter, and 
#' records the estimated parameter distributions.
#'
#' @param theta.off.diag.init A numeric vector representing the initial values 
#'   for the off-diagonal elements of the state transition matrix \code{Theta}.
#' @param theta.diag.init A numeric vector representing the initial values for 
#'   the diagonal elements of \code{Theta}.
#' @param ytn.sim A 3D array of dimensions \code{(N, p.dim, sim.num)} containing 
#'   the simulated observed data.
#' @param N An integer specifying the number of time points.
#' @param sim.num An integer specifying the number of simulation runs.
#' @param p.dim An integer specifying the number of observed variables.
#' @param m.dim An integer specifying the number of latent state variables.
#' @param tn.diff A numeric vector of length \code{N-1} representing the time 
#'   differences between consecutive observations.
#' @param theta.off.diag A numeric vector of the true off-diagonal elements of 
#'   \code{Theta} for comparison in confidence region calculations.
#' @param Z.par A numeric vector representing the true values of \code{Z} 
#'   for comparison in confidence region calculations.
#' @param H.diag A numeric vector representing the true diagonal elements of 
#'   \code{H} for comparison in confidence region calculations.
#' @param mu A numeric vector representing the true mean values for comparison.
#' @return A list containing:
#' \item{theta.matrix.dist}{An array storing the estimated \code{Theta} matrices across simulations.}
#' \item{C.matrix.dist}{An array storing the estimated state transition matrices \code{C}.}
#' \item{Z.matrix.dist}{An array storing the estimated observation matrices \code{Z}.}
#' \item{H.matrix.dist}{An array storing the estimated observation noise covariance matrices \code{H}.}
#' \item{mu.dist}{A matrix storing the estimated mean vectors.}
#' \item{Hessian.record}{An array storing the Hessian matrices from optimization.}
#' \item{CRegion}{A vector containing the confidence region values for overall parameter estimation.}
#' \item{CRegion.theta}{A vector containing the confidence region values for \code{Theta} parameters.}
#' \item{L.dist}{A vector storing the log-likelihood values from optimization.}
#' \item{result.para.dist}{A matrix storing the estimated parameter values.}
#' @details
#' The function follows these steps:
#' \itemize{
#'   \item Initializes the parameters based on the covariance of the data.
#'   \item Uses eigen decomposition to initialize the observation matrix \code{Z}.
#'   \item Runs the \code{optim()} function with the \code{L-BFGS-B} method to estimate parameters.
#'   \item Converts the estimated parameters into structured matrices using \code{para.to.model()}.
#'   \item Ensures sign consistency using \code{convert.Z.theta()}.
#'   \item Computes confidence regions based on the estimated Hessian matrices.
#' }
#' @examples
#' # Example usage (assuming required parameters and simulated data are predefined)
#' result <- estimation.OU(theta.off.diag.init, theta.diag.init, ytn.sim, N, sim.num, 
#'                         p.dim, m.dim, tn.diff, theta.off.diag, Z.par, H.diag, mu)
#' result$theta.matrix.dist  # Access the estimated Theta matrices
#' @export
estimation.OU <- function(theta.off.diag.init, theta.diag.init, ytn.sim, 
                          N, sim.num, p.dim, m.dim, tn.diff, 
                          theta.off.diag, Z.par, H.diag, mu) {
  # Initialize storage arrays
  C.matrix.dist <- array(NA, dim = c(m.dim, m.dim, sim.num))
  theta.matrix.dist <- array(NA, dim = c(m.dim, m.dim, sim.num))
  Z.matrix.dist <- array(NA, dim = c(p.dim, m.dim, sim.num))
  H.matrix.dist <- array(NA, dim = c(p.dim, p.dim, sim.num))
  L.dist <- rep(NA, sim.num)
  mu.dist <- matrix(NA, nrow = sim.num, ncol = p.dim)
  
  total.length.par <- p.dim * m.dim + p.dim + m.dim * (m.dim + 1) / 2 + p.dim
  result.para.dist <- matrix(NA, nrow = sim.num, ncol = total.length.par)
  Hessian.record <- array(NA, c(total.length.par, total.length.par, sim.num))
  
  CRegion <- rep(NA, sim.num)
  CRegion.theta <- rep(NA, sim.num)
  
  # Run estimation over multiple simulations
  for (i in 1:sim.num) {
    temp.sim <- matrix(ytn.sim[, , i], ncol = p.dim)
    
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
    theta.matrix.dist[, , i] <- convert.result$theta
    C.matrix.dist[, , i] <- convert.result$C.matrix
    Z.matrix.dist[, , i] <- convert.result$Z.matrix
    H.matrix.dist[, , i] <- result.mat$H.mat
    mu.dist[i, ] <- result.mat$mu
    Hessian.record[, , i] <- result.sim$hessian
    L.dist[i] <- result.sim$value
    result.para.dist[i, ] <- result.sim$par
    
    # Compute confidence region values
    param.diff <- result.sim$par - c(theta.off.diag, log(theta.diag), Z.par, log(H.diag), mu)
    CRegion[i] <- t(param.diff) %*% result.sim$hessian %*% param.diff
    
    param.diff.theta <- result.sim$par[1:(m.dim * (m.dim + 1) / 2)] - 
      c(theta.off.diag, log(theta.diag))
    CRegion.theta[i] <- t(param.diff.theta) %*% 
      result.sim$hessian[1:(m.dim * (m.dim + 1) / 2), 
                         1:(m.dim * (m.dim + 1) / 2)] %*% 
      param.diff.theta
  }
  
  # Return results
  return(list("theta.matrix.dist" = theta.matrix.dist, 
              "C.matrix.dist" = C.matrix.dist, 
              "Z.matrix.dist" = Z.matrix.dist, 
              "H.matrix.dist" = H.matrix.dist, 
              "mu.dist" = mu.dist, 
              "Hessian.record" = Hessian.record, 
              "CRegion" = CRegion, 
              "CRegion.theta" = CRegion.theta, 
              "L.dist" = L.dist, 
              "result.para.dist" = result.para.dist))
}
