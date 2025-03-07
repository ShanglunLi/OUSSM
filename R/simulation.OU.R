#' Simulate an Ornstein-Uhlenbeck Process
#'
#' This function simulates data from a multivariate Ornstein-Uhlenbeck (OU) process, 
#' generating both latent states and observed values over time.
#'
#' @param p.dim Integer. The number of observed variables.
#' @param m.dim Integer. The number of latent state variables.
#' @param mu Numeric vector of length \code{p.dim}. The mean of the observed process.
#' @param Z.matrix Numeric matrix of dimension \code{p.dim} x \code{m.dim}. The observation matrix.
#' @param H.matrix Numeric covariance matrix of dimension \code{p.dim} x \code{p.dim}. The observation noise covariance.
#' @param C.matrix List of transition matrices \eqn{C(t_n)} indexed by time differences.
#' @param Q.matrix List of process noise covariance matrices \eqn{Q(t_n)} indexed by time differences.
#' @param P1 Numeric covariance matrix of dimension \code{m.dim} x \code{m.dim}. The initial state covariance.
#' @param N Integer. The number of time points.
#' @param sim.num Integer. The number of simulations to run.
#' @param tn.diff Numeric vector of length \code{N-1}. The time differences between consecutive observations.
#' @return A list containing:
#'   \item{ytn}{An array of dimension \code{(N, p.dim, sim.num)}, representing the observed process.}
#'   \item{atn}{An array of dimension \code{(N, m.dim, sim.num)}, representing the latent process.}
#' @examples
#' # Example simulation (assuming required matrices are defined)
#' simulation.OU(p.dim = 2, m.dim = 2, mu = c(0, 0), 
#'               Z.matrix = matrix(c(1, 0, 0, 1), 2, 2), 
#'               H.matrix = diag(2), 
#'               C.matrix = list("1" = diag(2)), 
#'               Q.matrix = list("1" = diag(2)), 
#'               P1 = diag(2), 
#'               N = 10, sim.num = 5, 
#'               tn.diff = rep(1, 9))
#' @export
simulation.OU <- function(p.dim, m.dim, mu, Z.matrix, H.matrix, C.matrix, Q.matrix, P1, N, sim.num, tn.diff) {
  # Initialize arrays to store simulated values
  ytn <- array(NA, dim = c(N, p.dim, sim.num))
  atn <- array(NA, dim = c(N, m.dim, sim.num))
  
  # Loop over the number of simulations
  for (i in 1:sim.num) {
    # Generate observation noise
    epsilon.tn <- mvrnorm(n = N, mu = rep(0, p.dim), Sigma = H.matrix)
    
    # Initialize latent process noise
    eta.tn <- matrix(NA, nrow = N, ncol = m.dim)
    for (tn in 1:(N - 1)) {
      eta.tn[tn, ] <- mvrnorm(n = 1, mu = rep(0, m.dim), Sigma = Q.matrix[[as.character(tn.diff[tn])]])
    }
    
    # Generate initial state
    atn[1, , i] <- mvrnorm(n = 1, mu = rep(0, m.dim), Sigma = P1)
    
    # Simulate process over time
    for (tn in 1:N) {
      ytn[tn, , i] <- mu + Z.matrix %*% atn[tn, , i] + epsilon.tn[tn, ]
      if (tn < N) {
        atn[tn + 1, , i] <- C.matrix[[as.character(tn.diff[tn])]] %*% atn[tn, , i] + eta.tn[tn, ]
      }
    }
  }
  
  # Return results
  return(list("ytn" = ytn, "atn" = atn))
}
