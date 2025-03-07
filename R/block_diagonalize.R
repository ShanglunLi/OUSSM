#' Block Diagonalization of a Matrix
#'
#' This function performs block diagonalization of a given matrix \eqn{\Theta}, 
#' converting it into a real-valued block-diagonal form. It also applies the transformation 
#' to the matrices \eqn{Z} and \eqn{\Sigma}.
#'
#' @param Theta A square matrix to be block diagonalized.
#' @param Z A loading matrix associated with \eqn{\Theta}.
#' @param Sigma A covariance matrix associated with the process noise.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{block_matrix}: The block-diagonalized form of \eqn{\Theta}.
#'   \item \code{transformed_Z}: The transformed loading matrix.
#'   \item \code{transformed_Sigma}: The transformed covariance matrix.
#'   \item \code{transformation_matrix}: The real-valued transformation matrix used.
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Computes the eigenvalue decomposition of \eqn{\Theta}.
#'   \item Constructs a real-valued transformation matrix that block-diagonalizes \eqn{\Theta}.
#'   \item Transforms \eqn{Z} and \eqn{\Sigma} using this transformation.
#'   \item The block-diagonalized matrix consists of real values for real eigenvalues 
#'         and \eqn{2 \times 2} real blocks for complex-conjugate eigenvalues.
#' }
#'
#' @examples
#' # Example usage
#' Theta <- matrix(c(1, -2, 2, 3), nrow = 2)
#' Z <- matrix(c(0.5, 1.2, 0.8, -0.3), nrow = 2)
#' Sigma <- diag(2)
#' result <- blockDiagonalize(Theta, Z, Sigma)
#'
#' @export
block_diagonalize <- function(Theta, Z, Sigma) {
  library(Matrix)  # Load Matrix package for block diagonalization
  
  # Compute Eigen Decomposition
  eig <- eigen(Theta)
  values <- eig$values  # Eigenvalues
  vectors <- eig$vectors  # Eigenvectors (may be complex)
  
  # Construct the real block-diagonal transformation matrix
  used <- rep(FALSE, length(values))
  P_real <- matrix(0, nrow = nrow(Theta), ncol = ncol(Theta))
  blocks <- list()
  j <- 1  # Column index for P_real
  
  for (i in seq_along(values)) {
    if (!used[i]) {
      if (Im(values[i]) == 0) {
        # Real eigenvalue: Use the corresponding eigenvector
        P_real[, j] <- Re(vectors[, i])
        blocks <- c(blocks, list(matrix(Re(values[i]), 1, 1)))
        used[i] <- TRUE
        j <- j + 1
      } else {
        # Complex eigenvalue: Create a real-valued 2×2 block
        a <- Re(values[i])
        b <- Im(values[i])
        
        # Extract real and imaginary parts of the eigenvector
        v_real <- Re(vectors[, i])
        v_imag <- Im(vectors[, i])
        
        # Create a real-valued transformation
        P_real[, j] <- v_real
        P_real[, j + 1] <- v_imag
        
        # Construct the real 2x2 block
        complex_block <- matrix(c(a, b, -b, a), 2, 2)
        blocks <- c(blocks, list(complex_block))
        
        # Mark both complex-conjugate eigenvalues as used
        used[i] <- TRUE
        used[which(Im(values) == -b & Re(values) == a)] <- TRUE
        j <- j + 2
      }
    }
  }
  
  # Convert blocks into block-diagonal matrix
  block_diag_matrix <- as.matrix(bdiag(blocks))  # Convert to dense format
  
  # Compute the real transformation matrix P and its inverse
  P_inv <- solve(P_real)  # Inverse of P_real
  
  # Transform Z and Sigma
  Z_transformed <- Z %*% P_inv
  Sigma_transformed <- P_real %*% Sigma %*% t(P_real)
  
  return(list(
    block_matrix = block_diag_matrix,  # Block-diagonalized Theta
    transformed_Z = Z_transformed,     # Transformed Z
    transformed_Sigma = Sigma_transformed,  # Transformed Sigma
    transformation_matrix = P_real      # Real-valued transformation matrix
  ))
}
