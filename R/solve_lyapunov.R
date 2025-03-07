#' Solve the Lyapunov Equation
#'
#' This function solves the continuous Lyapunov equation \eqn{\Theta Q + Q \Theta^T = M}.
#'
#' @param Theta A square matrix.
#' @param M A square matrix of the same dimension as \code{Theta}.
#' @return A matrix \code{Q} that solves the equation.
#' @examples
#' Theta <- matrix(c(1, -2, 3, 4), nrow = 2)
#' M <- matrix(c(1, 0, 0, 1), nrow = 2)
#' solve_lyapunov(Theta, M)
#' @export
solve_lyapunov <- function(Theta, M) {
  d <- nrow(Theta)
  I_d <- diag(d)
  
  # Convert Lyapunov equation into linear system
  A_big <- kronecker(I_d, Theta) + kronecker(Theta, I_d)
  M_vec <- as.vector(M)
  
  # Solve for Q as a vector
  Q_vec <- solve(A_big, M_vec)
  
  # Reshape back into matrix form
  Q_mat <- matrix(Q_vec, nrow=d, ncol=d)
  return(Q_mat)
}