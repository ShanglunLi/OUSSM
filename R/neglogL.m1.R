#' Negative Log-Likelihood Function for 1D OUSSM
#'
#' This function computes the negative log-likelihood for a one-dimensional Ornstein-Uhlenbeck State-Space Model (OUSSM).
#'
#' @param para A numeric vector containing model parameters.
#' @param p.dim An integer specifying the number of observed variables.
#' @param N An integer specifying the total number of time points.
#' @param ytn A matrix of observed data with dimensions \code{(N, p.dim)}.
#' @param tn.diff A numeric vector representing the time gaps between observations.
#'
#' @return A numeric value representing the negative log-likelihood of the model.
#'
#' @details
#' The function follows these steps:
#' \itemize{
#'   \item Extracts and structures model parameters from \code{para}.
#'   \item Constructs the transition and noise matrices for the state-space model.
#'   \item Runs the Kalman filter to estimate the likelihood based on observed data.
#'   \item Computes the negative log-likelihood by summing the prediction errors over time.
#' }
#'
#' @examples
#' # Example usage (assuming required parameters and data are predefined)
#' result <- neglogL.m1(para, p.dim, N, ytn, tn.diff)
#'
#' @export
neglogL.m1 <- function(para, p.dim, N, ytn, tn.diff) {
  
  # Extract model parameters
  theta <- exp(para[1])
  para <- para[-c(1)]
  Z.mat <- matrix(para[1:p.dim], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:p.dim)]
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  mu <- para
  
  unique.t <- sort(unique(tn.diff))
  
  # Construct state transition matrices
  C.matrix.dic <- hash()
  for (t in unique.t) {
    C.matrix.dic[[as.character(t)]] <- exp(-theta * t)
  }
  
  # Construct noise covariance matrices
  H.matrix <- diag(H.diag)
  Q.matrix.dic <- hash()
  for (t in unique.t) {
    Q.matrix.dic[[as.character(t)]] <- (1 - exp(-2 * theta * t)) / (2 * theta)
  }
  
  # Initialize Kalman filter variables
  Ptn <- rep(NA, N)
  Ptn[1] <- 1 / (2 * theta)
  atn <- rep(0, N)
  Ftn <- array(NA, dim = c(p.dim, p.dim, N))
  Ftn[,,1] <- Z.mat %*% Ptn[1] %*% t(Z.mat) + H.matrix
  Ktn <- array(NA, dim = c(1, p.dim, (N - 1)))
  Ktn[,,1] <- C.matrix.dic[[as.character(1)]] %*% Ptn[1] %*% 
    t(Z.mat) %*% solve(Ftn[,,1])
  vtn <- matrix(NA, nrow = N, ncol = p.dim)
  vtn[1,] <- ytn[1,] - mu - Z.mat %*% atn[1]
  
  # Compute negative log-likelihood for the first observation
  if (p.dim > 1) {
    neg.logL <- (log(det(Ftn[,,1])) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
  } else {
    neg.logL <- (log(Ftn[,,1]) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
  }
  
  # Run Kalman filter iterations
  for (i in 1:(N - 1)) {
    atn[(i + 1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% 
      atn[i] + Ktn[,,i] %*% vtn[i,]
    Ptn[(i + 1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% 
      Ptn[i] %*% t(C.matrix.dic[[as.character(tn.diff[i])]] - 
                     Ktn[,,i] %*% Z.mat) + Q.matrix.dic[[as.character(tn.diff[i])]]
    Ftn[,,(i + 1)] <- Z.mat %*% Ptn[(i + 1)] %*% t(Z.mat) + H.matrix
    
    if (i < (N - 1)) {
      Ktn[,,(i + 1)] <- C.matrix.dic[[as.character(tn.diff[i + 1])]] %*% 
        Ptn[(i + 1)] %*% t(Z.mat) %*% solve(Ftn[,,(i + 1)])
    }
    
    vtn[(i + 1),] <- ytn[(i + 1),] - mu - Z.mat %*% atn[(i + 1)]
    
    if (p.dim > 1) {
      neg.logL <- neg.logL + (log(det(Ftn[,,(i + 1)])) + 
                                t(vtn[(i + 1),]) %*% solve(Ftn[,,(i + 1)]) %*% 
                                vtn[(i + 1),]) / 2
    } else {
      neg.logL <- neg.logL + (log(Ftn[,,(i + 1)]) + 
                                t(vtn[(i + 1),]) %*% solve(Ftn[,,(i + 1)]) %*% 
                                vtn[(i + 1),]) / 2
    }
  }
  
  return(neg.logL)
}
