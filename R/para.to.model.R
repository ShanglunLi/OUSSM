#' Convert Parameter Vector to Model Matrices
#'
#' This function converts a parameter vector into the corresponding model matrices for 
#' the Factor Ornstein-Uhlenbeck State-Space Model (OUSSM). It reconstructs the 
#' state transition matrix (\code{Theta}), observation matrix (\code{Z}), 
#' observation noise covariance (\code{H}), and mean vector (\code{mu}).
#'
#' @param para A numeric vector containing the model parameters, including off-diagonal 
#'   and diagonal elements of \code{Theta}, elements of \code{Z}, and diagonal elements 
#'   of \code{H}.
#' @param m.dim An integer specifying the number of latent state variables.
#' @param p.dim An integer specifying the number of observed variables.
#' @return A list containing:
#' \item{C.matrix}{State transition matrix, computed as \eqn{C = e^{-\Theta}}.}
#' \item{Z.mat}{Observation matrix of dimension \code{p.dim} x \code{m.dim}.}
#' \item{H.mat}{Observation noise covariance matrix (diagonal).}
#' \item{mu}{Mean vector of the observed process.}
#' \item{theta}{State transition matrix \code{Theta}, constructed with an antisymmetric 
#'   off-diagonal structure and decreasing diagonal elements.}
#' @details 
#' The function follows these steps:
#' \itemize{
#'   \item Extracts parameters from \code{para} and assigns them to respective matrices.
#'   \item Constructs \code{Theta} with off-diagonal elements forming an antisymmetric 
#'         structure and diagonal elements ensuring a decreasing order.
#'   \item Computes the state transition matrix \code{C} as \eqn{C = e^{-\Theta}}.
#'   \item Constructs the observation noise covariance \code{H} as a diagonal matrix.
#' }
#' @examples
#' # Example usage (assuming parameters are provided)
#' para <- c(runif(3), log(runif(2)), runif(6), log(runif(2)), runif(2))
#' model <- para.to.model(para, m.dim = 2, p.dim = 2)
#' model$C.matrix  # Access the state transition matrix
#' @export
para.to.model <- function(para, m.dim, p.dim) {
  # Extract off-diagonal elements of Theta
  theta.off.diag <- para[1:(m.dim * (m.dim - 1) / 2)]
  para <- para[-c(1:(m.dim * (m.dim - 1) / 2))]
  
  # Extract diagonal elements of Theta (exponentiated to ensure positivity)
  theta.diag <- exp(para[1:m.dim])
  para <- para[-c(1:m.dim)]
  
  # Extract observation matrix Z and reshape
  Z.mat <- matrix(para[1:(p.dim * m.dim)], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:(p.dim * m.dim))]
  
  # Extract and exponentiate observation noise diagonal elements
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  
  # Extract mean vector
  mu <- para
  
  # Construct Theta matrix
  theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
  theta.mat[upper.tri(theta.mat)] <- theta.off.diag
  theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
  
  # Ensure Theta's diagonal elements are decreasing (monotonicity constraint)
  theta.diag <- rev(cumsum(rev(theta.diag)))
  theta.mat <- theta.mat + diag(theta.diag)
  
  # Compute state transition matrix C
  C.matrix <- expm(-theta.mat)
  
  # Construct observation noise covariance matrix H
  H.matrix <- diag(H.diag)
  
  # Return a list of model components
  return(list("C.matrix" = C.matrix, 
              "Z.mat" = Z.mat, 
              "H.mat" = H.matrix, 
              "mu" = mu, 
              "theta" = theta.mat))
}