#' Compute the Negative Log-Likelihood for the Factor OUSSM Model
#'
#' This function computes the negative log-likelihood for the Factor Ornstein-Uhlenbeck 
#' State-Space Model (OUSSM) using the Kalman filter recursion.
#'
#' @param para A numeric vector of model parameters, including elements of the state 
#'   transition matrix (\code{Theta}), observation matrix (\code{Z}), noise covariance
#'   matrices (\code{H} and \code{Q}), and mean vector (\code{mu}).
#' @param p.dim An integer specifying the number of observed variables.
#' @param m.dim An integer specifying the number of latent state variables.
#' @param N An integer indicating the number of time points.
#' @param ytn A numeric matrix of observed data with dimensions \code{N x p.dim}.
#' @param tn.diff A numeric vector of length \code{N-1} representing time differences
#'   between consecutive observations.
#' @return A single numeric value representing the negative log-likelihood.
#' @details 
#' The function follows these steps:
#' \itemize{
#'   \item Extracts model parameters from \code{para} and reconstructs the matrices.
#'   \item Constructs the state transition matrix (\code{Theta}) ensuring it is a 
#'         sum of an antisymmetric matrix and a diagonal matrix with decreasing diagonal elements.
#'   \item Computes the state transition matrix exponentials (\code{C}) for unique time gaps.
#'   \item Computes process noise covariance matrices (\code{Q}) using the Lyapunov equation.
#'   \item Runs the Kalman filter to compute filtered state estimates and innovations.
#'   \item Accumulates the log-likelihood using the computed innovation covariances.
#' }
#' @examples
#' # Example usage (assuming required matrices and parameters are predefined)
#' para <- c(runif(3), log(runif(2)), runif(6), log(runif(2)), runif(2))
#' neglogL(para, p.dim = 2, m.dim = 2, N = 10, ytn = matrix(rnorm(20), 10, 2), tn.diff = rep(1, 9))
#' @export
neglogL <- function(para, p.dim, m.dim, N, ytn, tn.diff) {
  # Extract parameters for theta matrix (off-diagonal and diagonal)
  theta.off.diag <- para[1:(m.dim * (m.dim - 1) / 2)]
  para <- para[-c(1:(m.dim * (m.dim - 1) / 2))]
  
  # Diagonal elements of theta (exponentiated to ensure positivity)
  theta.diag <- exp(para[1:m.dim])
  para <- para[-c(1:m.dim)]
  
  # Observation matrix Z (reshaped)
  Z.mat <- matrix(para[1:(p.dim * m.dim)], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:(p.dim * m.dim))]
  
  # Observation noise covariance (diagonal elements, exponentiated)
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  
  # Mean vector
  mu <- para
  
  # Unique time differences
  unique.t <- sort(unique(tn.diff))
  
  # Construct theta matrix
  theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
  theta.mat[upper.tri(theta.mat)] <- theta.off.diag
  theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
  
  # Ensure decreasing diagonal elements
  theta.diag <- rev(cumsum(rev(theta.diag)))
  theta.mat <- theta.mat + diag(theta.diag)
  
  # Compute state transition matrices (C) for unique time differences
  C.matrix.dic <- hash()
  for (t in unique.t) {
    C.matrix.dic[[as.character(t)]] <- expm(-theta.mat * t)
  }
  
  # Construct observation noise covariance matrix (diagonal)
  H.matrix <- diag(H.diag)
  
  # Compute process noise covariance matrices (Q) using the Lyapunov solver
  Q.matrix.dic <- hash()
  for (t in unique.t) {
    M <- diag(m.dim) - C.matrix.dic[[as.character(t)]] %*% t(C.matrix.dic[[as.character(t)]])
    Q.matrix.dic[[as.character(t)]] <- solve_lyapunov(theta.mat, M)
  }
  
  # Initialize Kalman filter variables
  Ptn <- array(NA, dim = c(m.dim, m.dim, N))
  Ptn[, , 1] <- solve_lyapunov(theta.mat, diag(m.dim))  # Initial state covariance
  atn <- matrix(0, nrow = N, ncol = m.dim)  # Initial latent state
  
  # Initialize innovation covariance (F) and Kalman gain (K)
  Ftn <- array(NA, dim = c(p.dim, p.dim, N))
  Ftn[, , 1] <- Z.mat %*% Ptn[, , 1] %*% t(Z.mat) + H.matrix
  
  Ktn <- array(NA, dim = c(m.dim, p.dim, (N - 1)))
  Ktn[, , 1] <- C.matrix.dic[[as.character(1)]] %*% Ptn[, , 1] %*% t(Z.mat) %*% solve(Ftn[, , 1])
  
  # Compute initial innovation (v)
  vtn <- matrix(NA, nrow = N, ncol = p.dim)
  vtn[1, ] <- ytn[1, ] - mu - Z.mat %*% atn[1, ]
  
  # Compute initial negative log-likelihood
  if (p.dim > 1) {
    neg.logL <- (log(det(Ftn[, , 1])) + t(vtn[1, ]) %*% solve(Ftn[, , 1]) %*% vtn[1, ]) / 2
  } else {
    neg.logL <- (log(Ftn[, , 1]) + t(vtn[1, ]) %*% solve(Ftn[, , 1]) %*% vtn[1, ]) / 2
  }
  
  # Run Kalman filter recursion
  for (i in 1:(N - 1)) {
    atn[(i + 1), ] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% atn[i, ] + Ktn[, , i] %*% vtn[i, ]
    Ptn[, , (i + 1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% 
      Ptn[, , i] %*% 
      t(C.matrix.dic[[as.character(tn.diff[i])]] - Ktn[, , i] %*% Z.mat) + 
      Q.matrix.dic[[as.character(tn.diff[i])]]
    
    Ftn[, , (i + 1)] <- Z.mat %*% Ptn[, , (i + 1)] %*% t(Z.mat) + H.matrix
    
    if (i < (N - 1)) {
      Ktn[, , (i + 1)] <- C.matrix.dic[[as.character(tn.diff[i + 1])]] %*% 
        Ptn[, , (i + 1)] %*% 
        t(Z.mat) %*% solve(Ftn[, , (i + 1)])
    }
    
    vtn[(i + 1), ] <- ytn[(i + 1), ] - mu - Z.mat %*% atn[(i + 1), ]
    
    # Update negative log-likelihood
    if (p.dim > 1) {
      neg.logL <- neg.logL + (log(det(Ftn[, , (i + 1)])) + 
                                t(vtn[(i + 1), ]) %*% solve(Ftn[, , (i + 1)]) %*% vtn[(i + 1), ]) / 2
    } else {
      neg.logL <- neg.logL + (log(Ftn[, , (i + 1)]) + 
                                t(vtn[(i + 1), ]) %*% solve(Ftn[, , (i + 1)]) %*% vtn[(i + 1), ]) / 2
    }
  }
  
  return(neg.logL)
}

