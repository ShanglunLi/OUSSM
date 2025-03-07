
create_model <- function(theta.off.diag, theta.diag,  Z.par, H.diag, tn.diff, m.dim){
  p.dim <- length(Z.par) / m.dim
  
  unique.t <- sort(unique(tn.diff))
  Z.matrix <- matrix(Z.par, nrow = p.dim, byrow = TRUE)
  
  theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
  theta.mat[upper.tri(theta.mat)] <- theta.off.diag
  theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
  theta.mat <- theta.mat + diag(theta.diag)
  
  C.matrix.dic <- hash()
  for (t in unique.t) {
    C.matrix <- expm(-theta.mat*t)
    C.matrix.dic[[as.character(t)]] <- C.matrix
  }
  
  H.matrix <- diag(H.diag) 
  
  Q.matrix.dic <- hash()
  for (t in unique.t) {
    M <- diag(m.dim) - C.matrix.dic[[as.character(t)]] %*% t(C.matrix.dic[[as.character(t)]])
    Q.matrix.dic[[as.character(t)]] <- solve_lyapunov(theta.mat, M)
  }
  
  P1 <- solve_lyapunov(theta.mat, diag(m.dim))
  
  result <- list("Z.mat" = Z.matrix, "Q.mat" =  Q.matrix.dic, "theta" = theta.mat, "C.matrix" = C.matrix.dic, "H.mat" = H.matrix, "P1.mat" = P1)
  return(result)
}