#' Create a State-Space Model
#'
#' This function constructs a state-space model by computing various matrices, including 
#' the transition matrix (\eqn{\Theta}), observation matrix (\eqn{Z}), and noise covariance matrices.
#'
#' @param theta.off.diag A vector of off-diagonal elements for the \code{Theta} matrix (upper triangular part).
#' @param theta.diag A vector of diagonal elements for the \code{Theta} matrix.
#' @param Z.par A vector specifying the elements of the observation matrix \code{Z}.
#' @param H.diag A vector specifying the diagonal elements of the observation noise covariance matrix \code{H}.
#' @param tn.diff A numeric vector of time differences for matrix exponentiation.
#' @param m.dim An integer specifying the number of latent states.
#' @return A list containing:
#' \itemize{
#'   \item \code{Z.mat} - The observation matrix.
#'   \item \code{Q.mat} - A dictionary of process noise covariance matrices.
#'   \item \code{theta} - The state transition matrix.
#'   \item \code{C.matrix} - A dictionary of transition matrices for different time steps.
#'   \item \code{H.mat} - The observation noise covariance matrix.
#'   \item \code{P1.mat} - The stationary diffusion error covariance matrix.
#' }
#' @examples
#' theta.off.diag <- c(0.5, -0.3, 0.2)
#' theta.diag <- c(1.0, 0.8, 0.6)
#' Z.par <- c(1, 0.5, 0.3, 0.4, 1, 0.2)
#' H.diag <- c(0.1, 0.1)
#' tn.diff <- c(1, 2, 3)
#' m.dim <- 3
#' create_model(theta.off.diag, theta.diag, Z.par, H.diag, tn.diff, m.dim)
#' @export
create_model <- function(theta.off.diag, theta.diag, Z.par, H.diag, tn.diff, m.dim) {
  library(expm)  # For matrix exponentiation
  library(hash)  # For storing matrices in a dictionary
  
  p.dim <- length(Z.par) / m.dim
  
  # Create the Z matrix
  Z.matrix <- matrix(Z.par, nrow = p.dim, byrow = TRUE)
  
  # Construct the Theta matrix
  theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
  theta.mat[upper.tri(theta.mat)] <- theta.off.diag
  theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
  theta.mat <- theta.mat + diag(theta.diag)
  
  # Create dictionaries for C and Q matrices
  C.matrix.dic <- hash()
  Q.matrix.dic <- hash()
  
  unique.t <- sort(unique(tn.diff))
  
  # Compute matrix exponentials and store in dictionary
  for (t in unique.t) {
    C.matrix <- expm(-theta.mat * t)
    C.matrix.dic[[as.character(t)]] <- C.matrix
  }
  
  # Construct H matrix
  H.matrix <- diag(H.diag)
  
  # Compute Q matrices using Lyapunov equation
  for (t in unique.t) {
    M <- diag(m.dim) - C.matrix.dic[[as.character(t)]] %*% t(C.matrix.dic[[as.character(t)]])
    Q.matrix.dic[[as.character(t)]] <- solve_lyapunov(theta.mat, M)
  }
  
  # Compute P1 matrix
  P1 <- solve_lyapunov(theta.mat, diag(m.dim))
  
  # Return a list of all computed matrices
  result <- list(
    "Z.mat" = Z.matrix, 
    "Q.mat" = Q.matrix.dic, 
    "theta" = theta.mat, 
    "C.matrix" = C.matrix.dic, 
    "H.mat" = H.matrix, 
    "P1.mat" = P1
  )
  
  return(result)
}
