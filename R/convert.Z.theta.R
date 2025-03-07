#' Convert Z and Theta Matrices for Sign Consistency
#'
#' This function modifies the \code{Z} and \code{Theta} matrices to ensure sign consistency
#' in the Factor Ornstein-Uhlenbeck State-Space Model (OUSSM). It constructs a transformation 
#' matrix \eqn{M} that normalizes the sign of the first row of \code{Z}, then applies it to 
#' both \code{Theta} and \code{Z}. The state transition matrix \code{C} is also updated.
#'
#' @param theta.mat A square matrix representing the state transition matrix \code{Theta}.
#' @param Z.matrix A matrix representing the observation loading matrix \code{Z}.
#' @return A list containing:
#' \item{theta.mat}{Updated state transition matrix with sign consistency applied.}
#' \item{Z.matrix}{Updated observation matrix with consistent sign conventions.}
#' \item{C.matrix}{State transition matrix computed as \eqn{C = e^{-\Theta}}.}
#' @details
#' The function applies the following transformation:
#' \itemize{
#'   \item Constructs a diagonal matrix \code{M} that normalizes the first row of \code{Z}.
#'   \item Applies \code{M} to both \code{Z} and \code{Theta} to maintain sign consistency.
#'   \item Computes the updated state transition matrix \code{C = expm(-Theta)}.
#' }
#' @examples
#' # Example usage
#' theta.mat <- matrix(c(1, -0.5, 0.5, 2), nrow = 2)
#' Z.matrix <- matrix(c(0.8, -1.2, 1.5, -0.7), nrow = 2)
#' converted <- convert.Z.theta(theta.mat, Z.matrix)
#' converted$theta.mat  # Access transformed Theta matrix
#' @export
convert.Z.theta <- function(theta.mat, Z.matrix) {
  # Construct transformation matrix M based on the first row of Z
  M.mat <- diag(Z.matrix[1,] / abs(Z.matrix[1,]))  # Ensures sign consistency
  
  # Apply transformation to Z and Theta
  Z.matrix.new <- Z.matrix %*% M.mat
  theta.mat.new <- M.mat %*% theta.mat %*% M.mat
  
  # Compute updated state transition matrix C
  C.matrix <- expm(-theta.mat.new)
  
  # Return transformed matrices
  return(list("theta.mat" = theta.mat.new, 
              "Z.matrix" = Z.matrix.new, 
              "C.matrix" = C.matrix))
}
