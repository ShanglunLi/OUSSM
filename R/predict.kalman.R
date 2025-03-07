#' Kalman Filter-Based Prediction for Factor OUSSM
#'
#' This function applies the Kalman filter to predict future observations 
#' using a trained Factor Ornstein-Uhlenbeck State-Space Model (OUSSM).
#'
#' @param para A numeric vector of estimated model parameters.
#' @param p.dim An integer specifying the number of observed variables.
#' @param m.dim An integer specifying the number of latent state variables.
#' @param N An integer specifying the total number of time points.
#' @param ytn A matrix of observed data with dimensions \code{(N, p.dim)}.
#' @param tn.diff A numeric vector representing the time gaps between observations.
#' @param N.test An integer specifying the number of time points for testing.
#'
#' @return A list containing:
#' \item{ytn.pred}{A matrix of predicted values for the test set.}
#' \item{ytn.upper}{A matrix representing the upper bounds of the 95\% confidence interval.}
#' \item{ytn.lower}{A matrix representing the lower bounds of the 95\% confidence interval.}
#'
#' @details
#' The function follows these steps:
#' \itemize{
#'   \item Extracts and structures model parameters from \code{para}.
#'   \item Constructs state transition and noise matrices.
#'   \item Runs the Kalman filter to estimate latent states and one-step predictions.
#'   \item Computes confidence intervals for the predictions.
#' }
#'
#' @examples
#' # Example usage (assuming required parameters and data are predefined)
#' result <- predict.kalman(para, p.dim, m.dim, N, ytn, tn.diff, N.test)
#' result$ytn.pred  # Access predicted values
#'
#' @export
predict.kalman <- function(para, p.dim, m.dim, N, ytn, tn.diff, N.test) {
  
  # Extract model parameters
  theta.off.diag <- para[1:(m.dim * (m.dim - 1) / 2)]
  para <- para[-c(1:(m.dim * (m.dim - 1) / 2))]
  theta.diag <- exp(para[1:m.dim])
  para <- para[-c(1:m.dim)]
  Z.mat <- matrix(para[1:(p.dim * m.dim)], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:(p.dim * m.dim))]
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  mu <- para
  
  # Define number of training points
  N.train <- N - N.test
  
  # Unique time differences
  unique.t <- sort(unique(tn.diff))
  
  # Construct Theta matrix
  theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
  theta.mat[upper.tri(theta.mat)] <- theta.off.diag
  theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
  theta.diag <- rev(cumsum(rev(theta.diag)))
  theta.mat <- theta.mat + diag(theta.diag)
  
  # Construct state transition matrices
  C.matrix.dic <- hash()
  for (t in unique.t) {
    C.matrix <- expm(-theta.mat * t)
    C.matrix.dic[[as.character(t)]] <- C.matrix
  }
  
  # Construct covariance matrices
  H.matrix <- diag(H.diag)
  Q.matrix.dic <- hash()
  for (t in unique.t) {
    M <- diag(m.dim) - C.matrix.dic[[as.character(t)]] %*% 
      t(C.matrix.dic[[as.character(t)]])
    Q.matrix.dic[[as.character(t)]] <- solve_lyapunov(theta.mat, M)
  }
  
  # Initialize Kalman filter variables
  Ptn <- array(NA, dim = c(m.dim, m.dim, N))
  Ptn[,,1] <- solve_lyapunov(theta.mat, diag(m.dim))
  atn <- matrix(0, nrow = N, ncol = m.dim)
  Ftn <- array(NA, dim = c(p.dim, p.dim, N))
  Ftn[,,1] <- Z.mat %*% Ptn[,,1] %*% t(Z.mat) + H.matrix
  Ktn <- array(NA, dim = c(m.dim, p.dim, (N - 1)))
  Ktn[,,1] <- C.matrix.dic[[as.character(1)]] %*% Ptn[,,1] %*% 
    t(Z.mat) %*% solve(Ftn[,,1])
  vtn <- matrix(NA, nrow = N, ncol = p.dim)
  vtn[1,] <- ytn[1,] - mu - Z.mat %*% atn[1,]
  
  # Initialize prediction storage
  ytn.pred <- matrix(NA, nrow = N.test, ncol = p.dim)
  ytn.upper <- matrix(NA, nrow = N.test, ncol = p.dim)
  ytn.lower <- matrix(NA, nrow = N.test, ncol = p.dim)
  
  # Run Kalman filter
  for (i in 1:(N - 1)) {
    atn[(i + 1),] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% 
      atn[i,] + Ktn[,,i] %*% vtn[i,]
    Ptn[,,(i + 1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% 
      Ptn[,,i] %*% t(C.matrix.dic[[as.character(tn.diff[i])]] - 
                       Ktn[,,i] %*% Z.mat) + Q.matrix.dic[[as.character(tn.diff[i])]]
    Ftn[,,(i + 1)] <- Z.mat %*% Ptn[,,(i + 1)] %*% t(Z.mat) + H.matrix
    
    # Generate predictions for test set
    if (i >= N.train) {
      ytn.pred[(i - N.train + 1),] <- mu + Z.mat %*% atn[(i + 1),]
      for (F.i in 1:p.dim) {
        ytn.upper[(i - N.train + 1), F.i] <- 
          ytn.pred[(i - N.train + 1), F.i] + 
          1.96 * sqrt(Ftn[F.i, F.i, (i - N.train + 1)])
        ytn.lower[(i - N.train + 1), F.i] <- 
          ytn.pred[(i - N.train + 1), F.i] - 
          1.96 * sqrt(Ftn[F.i, F.i, (i - N.train + 1)])
      }
    }
    
    if (i < (N - 1)) {
      Ktn[,,(i + 1)] <- C.matrix.dic[[as.character(tn.diff[i + 1])]] %*% 
        Ptn[,,(i + 1)] %*% t(Z.mat) %*% solve(Ftn[,,(i + 1)])
    }
    
    vtn[(i + 1),] <- ytn[(i + 1),] - mu - Z.mat %*% atn[(i + 1),]
  }
  
  # Return predictions and confidence intervals
  return(list("ytn.pred" = ytn.pred, 
              "ytn.upper" = ytn.upper, 
              "ytn.lower" = ytn.lower))
}
