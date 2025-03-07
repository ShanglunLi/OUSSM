#' Convert Parameter Vector to Model Matrices (1D OUSSM)
#'
#' This function transforms a vector of model parameters into structured model matrices 
#' for a one-dimensional Ornstein-Uhlenbeck State-Space Model (OUSSM).
#'
#' @param para A numeric vector containing model parameters.
#' @param p.dim An integer specifying the number of observed variables.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{C.matrix}: The state transition matrix.
#'   \item \code{Z.mat}: The loading matrix.
#'   \item \code{H.mat}: The measurement noise covariance matrix.
#'   \item \code{mu}: The mean vector of observations.
#'   \item \code{theta}: The estimated mean-reverting rate.
#' }
#'
#' @details
#' The function follows these steps:
#' \itemize{
#'   \item Extracts and applies constraints (exponential transformation) on \code{theta} and \code{H.diag}.
#'   \item Constructs the transition matrix \code{C.matrix} using the estimated mean-reverting rate.
#'   \item Reshapes and assigns the extracted parameter values to model components.
#' }
#'
#' @examples
#' # Example usage (assuming required parameters are predefined)
#' para <- c(log(0.5), 0.3, 0.2, log(0.1), log(0.05), 0.5, 0.4)
#' p.dim <- 2
#' model <- para.to.model.m1(para, p.dim)
#'
#' @export
para.to.model.m1 <- function(para, p.dim) {
  
  # Extract mean-reverting rate and apply transformation
  theta <- exp(para[1])
  para <- para[-c(1)]
  
  # Construct loading matrix
  Z.mat <- matrix(para[1:p.dim], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:p.dim)]
  
  # Extract and apply transformation on measurement noise covariance
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  
  # Extract mean vector
  mu <- para
  
  # Construct state transition matrix
  C.matrix <- exp(-theta)
  
  # Construct measurement noise covariance matrix
  H.matrix <- diag(H.diag)
  
  # Return structured model components
  result <- list("C.matrix" = C.matrix, 
                 "Z.mat" = Z.mat, 
                 "H.mat" = H.matrix, 
                 "mu" = mu, 
                 "theta" = theta)
  
  return(result)
}
