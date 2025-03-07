##### OU process simulation #####
# create the matrices of the model


create_model_small_gap <- function(theta.off.diag, theta.diag,  Z.par, H.diag, tn.diff, m.dim){
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

# Description: this function is used to simulate OU process with p,m dimension
# model:
# ytn = mu + A %*% atn + epsilontn, epsilontn ~ N(0, H)
# atn+1 = Ctn %*% atn + etatn, eta ~ N(0, Qtn)
# p.dim is the length of ytn, m.dim is the length of atn
simulation.OU <- function(p.dim, m.dim, mu, Z.matrix, H.matrix, C.matrix, Q.matrix, P1, N, sim.num, tn.diff) {
  ytn <- array(NA, dim = c(N,p.dim,sim.num))
  atn <- array(NA, dim = c(N,m.dim,sim.num))
  
  # simulate ytn
  for (i in 1:sim.num) {
    # simulate epsilon_tn
    epsilon.tn <- mvrnorm(n = N, mu = rep(0, p.dim), Sigma = H.matrix)
    
    # simulate eta_tn
    eta.tn <- matrix(NA, nrow = N, ncol = m.dim)
    for (tn in 1:(N-1)) {
      
      eta.tn[tn,] <- mvrnorm(n = 1, mu = rep(0, m.dim), Sigma = Q.matrix[[as.character(tn.diff[tn])]])
    }
    
    atn[1,,i] <- mvrnorm(n = 1, mu = rep(0, m.dim), Sigma = P1) # at0
    for (tn in 1:N) {
      ytn[tn,,i] <- mu + Z.matrix %*% atn[tn,,i] + epsilon.tn[tn,]
      if (tn < N) {
        atn[(tn + 1),,i] <- C.matrix[[as.character(tn.diff[tn])]] %*% atn[tn,,i] + eta.tn[tn,]
      }
    }
  }
  result <- list("ytn" = ytn, "atn" = atn)
  
  return(result)
}

### own neglogL
neglogL <- function(para, p.dim, m.dim, N, ytn, tn.diff) {
  theta.off.diag <- para[1:(m.dim*(m.dim-1)/2)]
  para <- para[-c(1:(m.dim*(m.dim-1)/2))]
  theta.diag <- exp(para[1:m.dim])
  para <- para[-c(1:m.dim)]
  Z.mat <- matrix(para[1:(p.dim*m.dim)], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:(p.dim*m.dim))]
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  mu <- para

  unique.t <- sort(unique(tn.diff))
  theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
  theta.mat[upper.tri(theta.mat)] <- theta.off.diag
  theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
  theta.diag <- rev(cumsum(rev(theta.diag)))
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
  
  Ptn <- array(NA, dim = c(m.dim,m.dim,N))
  Ptn[,,1] <- solve_lyapunov(theta.mat, diag(m.dim))
  atn <- matrix(0, nrow = N, ncol = m.dim)
  
  Ftn <- array(NA, dim = c(p.dim,p.dim,N))
  Ftn[,,1] <- Z.mat %*% Ptn[,,1] %*% t(Z.mat) + H.matrix

  Ktn <- array(NA, dim = c(m.dim,p.dim,(N-1)))
  Ktn[,,1] <- C.matrix.dic[[as.character(1)]] %*% Ptn[,,1] %*% t(Z.mat) %*% solve(Ftn[,,1])
  
  vtn <- matrix(NA, nrow = N, ncol = p.dim)
  vtn[1,] <- ytn[1,] - mu - Z.mat %*% atn[1,]

  # !!!!!! -(N + 1)*p/2 * log(2 * pi)
  if (p.dim > 1) {
    neg.logL <- (log(det(Ftn[,,1])) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
  }
  else {
    neg.logL <- (log(Ftn[,,1]) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
  }
  
  for (i in 1:(N-1)) {
    atn[(i+1),] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% atn[i,] + Ktn[,,i] %*% vtn[i,]
    Ptn[,,(i+1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% Ptn[,,i] %*% t(C.matrix.dic[[as.character(tn.diff[i])]] - Ktn[,,i] %*% Z.mat) + Q.matrix.dic[[as.character(tn.diff[i])]]
    Ftn[,,(i+1)] <- Z.mat %*% Ptn[,,(i+1)] %*% t(Z.mat) + H.matrix
    if (i < (N-1)) {
      Ktn[,,(i+1)] <- C.matrix.dic[[as.character(tn.diff[i+1])]] %*% Ptn[,,(i+1)] %*% t(Z.mat) %*% solve(Ftn[,,(i+1)])
    }
    vtn[(i+1),] <- ytn[(i+1),] - mu - Z.mat %*% atn[(i+1),]
    if (p.dim > 1) {
      neg.logL <- neg.logL + (log(det(Ftn[,,(i+1)])) + t(vtn[(i+1),]) %*% solve(Ftn[,,(i+1)]) %*% vtn[(i+1),]) / 2
    }
    else {
      neg.logL <- neg.logL + (log(Ftn[,,(i+1)]) + t(vtn[(i+1),]) %*% solve(Ftn[,,(i+1)]) %*% vtn[(i+1),]) / 2
    }
    
  }
  
  return(neg.logL)
}

para.to.model <- function(para, m.dim, p.dim){
  theta.off.diag <- para[1:(m.dim*(m.dim-1)/2)]
  para <- para[-c(1:(m.dim*(m.dim-1)/2))]
  theta.diag <- exp(para[1:m.dim])
  para <- para[-c(1:m.dim)]
  Z.mat <- matrix(para[1:(p.dim*m.dim)], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:(p.dim*m.dim))]
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  mu <- para
  
  theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
  theta.mat[upper.tri(theta.mat)] <- theta.off.diag
  theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
  theta.diag <- rev(cumsum(rev(theta.diag)))
  theta.mat <- theta.mat + diag(theta.diag)
  
  C.matrix <- expm(-theta.mat)
  
  H.matrix <- diag(H.diag) 
  
  result <- list("C.matrix" = C.matrix, "Z.mat" = Z.mat, "H.mat" = H.matrix, "mu" = mu, "theta" = theta.mat)
  return(result)
}

estimation.OU <- function(theta.off.diag.init, theta.diag.init, ytn.sim, N, sim.num, p.dim, m.dim, tn.diff, theta.off.diag, Z.par, H.diag, mu) {
  C.matrix.dist <- array(NA, dim = c(m.dim, m.dim, sim.num))
  theta.matrix.dist <- array(NA, dim = c(m.dim, m.dim, sim.num))
  Z.matrix.dist <- array(NA, dim = c(p.dim, m.dim, sim.num))
  H.matrix.dist <- array(NA, dim = c(p.dim, p.dim, sim.num))
  L.dist <- rep(NA, sim.num)
  mu.dist <- matrix(NA, nrow = sim.num, ncol = p.dim)
  total.length.par <- p.dim*m.dim + p.dim + m.dim*(m.dim+1)/2 + p.dim
  result.para.dist <- matrix(NA, nrow = sim.num, ncol = total.length.par)
  Hessian.record <- array(NA, c(total.length.par, total.length.par, sim.num))
  CRegion <- rep(NA, sim.num)
  CRegion.theta <- rep(NA, sim.num)
  
  for (i in 1:sim.num) {
    temp.sim <- matrix(ytn.sim[,,i], ncol = p.dim)
    
    M <- cov(temp.sim[1:(N-1),], temp.sim[2:N,])
    M <- eigen(M + t(M))
    if (p.dim >= m.dim) {
      Z.init <- M$vectors[,1:m.dim]
    }
    if (p.dim < m.dim) {
      Z.init <- cbind(M$vectors, matrix(1,nrow = p.dim, ncol = (m.dim-p.dim)))
    }
    Z.init <- c(t(Z.init))
    H.par.init <- diag(cov(temp.sim)) / 3
    mu.init <- colMeans(temp.sim)
    init.par <- c(theta.off.diag.init, log(theta.diag.init), Z.init, H.par.init, mu.init)
    
    result.sim <- optim(init.par, p.dim = p.dim, m.dim = m.dim, N = N, ytn = temp.sim, tn.diff = tn.diff, neglogL, method = "L-BFGS-B", hessian = TRUE)
    result.mat <- para.to.model(result.sim$par, m.dim, p.dim)
    
    convert.result <- convert.Z.theta(result.mat$theta, result.mat$Z.mat)
    theta.matrix.dist[,,i] <- convert.result$theta
    C.matrix.dist[,,i] <- convert.result$C.matrix
    Z.matrix.dist[,,i] <- convert.result$Z.matrix
    H.matrix.dist[,,i] <- result.mat$H.mat
    mu.dist[i,] <- result.mat$mu
    Hessian.record[,,i] <- result.sim$hessian
    L.dist[i] <- result.sim$value
    result.para.dist[i,] <- result.sim$par

    # Calculate Confidence region
    CRegion[i] <- t(result.sim$par - c(theta.off.diag, log(theta.diag), Z.par, log(H.diag), mu)) %*% (result.sim$hessian) %*% (result.sim$par - c(theta.off.diag, log(theta.diag), Z.par, log(H.diag), mu))
    CRegion.theta[i] <- t(result.sim$par[1:(m.dim*(m.dim+1)/2)] - c(theta.off.diag, log(theta.diag))) %*% (result.sim$hessian[1:(m.dim*(m.dim+1)/2),1:(m.dim*(m.dim+1)/2)]) %*% (result.sim$par[1:(m.dim*(m.dim+1)/2)] - c(theta.off.diag, log(theta.diag)))
  }
  
  result <- list("theta.matrix.dist" = theta.matrix.dist, "C.matrix.dist" = C.matrix.dist, "Z.matrix.dist" = Z.matrix.dist, "H.matrix.dist" = H.matrix.dist, "mu.dist" = mu.dist, "Hessian.record" = Hessian.record, "CRegion" = CRegion, "CRegion.theta" = CRegion.theta, "L.dist" = L.dist, "result.para.dist" = result.para.dist)
  
  return(result)
}

convert.Z.theta <- function(theta.mat, Z.matrix) {
  M.mat <- diag(Z.matrix[1,] / abs(Z.matrix[1,]))
  Z.matrix.new <- Z.matrix %*% M.mat
  theta.mat.new <- M.mat %*% theta.mat %*% M.mat
  C.matrix <- expm(-theta.mat.new)
  
  result <- list("theta.mat" = theta.mat.new, "Z.matrix" = Z.matrix.new, "C.matrix" = C.matrix)
  
  return(result)
}

estimation.OU.single <- function(theta.off.diag.init, theta.diag.init, ytn.sim, N, p.dim, m.dim, tn.diff) {
  temp.sim <- as.matrix(ytn.sim)
  
  M <- cov(temp.sim[1:(N-1),], temp.sim[2:N,])
  M <- eigen(M + t(M))
  if (p.dim >= m.dim) {
    Z.init <- M$vectors[,1:m.dim]
  }
  if (p.dim < m.dim) {
    Z.init <- cbind(M$vectors, matrix(1,nrow = p.dim, ncol = (m.dim-p.dim)))
  }
  Z.init <- c(t(Z.init))
  H.par.init <- diag(cov(temp.sim)) / 3
  mu.init <- colMeans(temp.sim)
  init.par <- c(theta.off.diag.init, log(theta.diag.init), Z.init, H.par.init, mu.init)
  
  result.sim <- optim(init.par, p.dim = p.dim, m.dim = m.dim, N = N, ytn = temp.sim, tn.diff = tn.diff, neglogL, method = "L-BFGS-B", hessian = TRUE)
  result.mat <- para.to.model(result.sim$par, m.dim, p.dim)
  
  convert.result <- convert.Z.theta(result.mat$theta, result.mat$Z.mat)
  theta.matrix.dist <- convert.result$theta
  C.matrix.dist <- convert.result$C.matrix
  Z.matrix.dist <- convert.result$Z.matrix
  H.matrix.dist <- result.mat$H.mat
  mu.dist <- result.mat$mu
  Hessian.record <- result.sim$hessian
  L.dist <- result.sim$value
  result.para.dist <- result.sim$par
  
  AIC.value <- 2 * L.dist + 2 * length(result.para.dist)
  BIC.value <- 2 * L.dist + log(N) * length(result.para.dist)
  
  result <- list("theta.matrix.dist" = theta.matrix.dist, "C.matrix.dist" = C.matrix.dist, "Z.matrix.dist" = Z.matrix.dist, "H.matrix.dist" = H.matrix.dist, "mu.dist" = mu.dist, "Hessian.record" = Hessian.record, 
                 "L.dist" = L.dist, "result.para.dist" = result.para.dist, "AIC" = AIC.value, "BIC" = BIC.value)
  
  return(result)
}

predict.kalman <- function(para, p.dim, m.dim, N, ytn, tn.diff, N.test) {
  theta.off.diag <- para[1:(m.dim*(m.dim-1)/2)]
  para <- para[-c(1:(m.dim*(m.dim-1)/2))]
  theta.diag <- exp(para[1:m.dim])
  para <- para[-c(1:m.dim)]
  Z.mat <- matrix(para[1:(p.dim*m.dim)], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:(p.dim*m.dim))]
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  mu <- para
  
  N.train <- N - N.test

  unique.t <- sort(unique(tn.diff))
  theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
  theta.mat[upper.tri(theta.mat)] <- theta.off.diag
  theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
  theta.diag <- rev(cumsum(rev(theta.diag)))
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
  
  Ptn <- array(NA, dim = c(m.dim,m.dim,N))
  Ptn[,,1] <- solve_lyapunov(theta.mat, diag(m.dim))
  atn <- matrix(0, nrow = N, ncol = m.dim)
  
  Ftn <- array(NA, dim = c(p.dim,p.dim, N))
  Ftn[,,1] <- Z.mat %*% Ptn[,,1] %*% t(Z.mat) + H.matrix
  
  Ktn <- array(NA, dim = c(m.dim,p.dim,(N -1)))
  Ktn[,,1] <- C.matrix.dic[[as.character(1)]] %*% Ptn[,,1] %*% t(Z.mat) %*% solve(Ftn[,,1])
  
  vtn <- matrix(NA, nrow = N, ncol = p.dim)
  vtn[1,] <- ytn[1,] - mu - Z.mat %*% atn[1,]
  
  ytn.pred <- matrix(NA, nrow = N.test, ncol = p.dim)
  ytn.upper <- matrix(NA, nrow = N.test, ncol = p.dim)
  ytn.lower <- matrix(NA, nrow = N.test, ncol = p.dim)
  for (i in 1:(N -1)) {
    atn[(i+1),] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% atn[i,] + Ktn[,,i] %*% vtn[i,]
    Ptn[,,(i+1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% Ptn[,,i] %*% t(C.matrix.dic[[as.character(tn.diff[i])]] - Ktn[,,i] %*% Z.mat) + Q.matrix.dic[[as.character(tn.diff[i])]]
    Ftn[,,(i+1)] <- Z.mat %*% Ptn[,,(i+1)] %*% t(Z.mat) + H.matrix
    if (i >= N.train) {
      ytn.pred[(i-N.train+1),] <- mu + Z.mat %*% atn[(i+1),]
      for (F.i in 1:p.dim) {
        ytn.upper[(i-N.train+1),F.i] <- ytn.pred[(i-N.train+1),F.i] + 1.96 * sqrt(Ftn[F.i,F.i,(i-N.train+1)])
        ytn.lower[(i-N.train+1),F.i] <- ytn.pred[(i-N.train+1),F.i] - 1.96 * sqrt(Ftn[F.i,F.i,(i-N.train+1)])
      }
    }
    if (i < (N -1)) {
      Ktn[,,(i+1)] <- C.matrix.dic[[as.character(tn.diff[i+1])]] %*% Ptn[,,(i+1)] %*% t(Z.mat) %*% solve(Ftn[,,(i+1)])
    }
    vtn[(i+1),] <- ytn[(i+1),] - mu - Z.mat %*% atn[(i+1),]
    
  }
  
  result <- list("ytn.pred" = ytn.pred, "ytn.upper" = ytn.upper, "ytn.lower" = ytn.lower)
  
  return(result)
}

simulation.trace <- function(p.dim, m.dim, mu, Z.matrix, H.matrix, C.matrix, Q.matrix, P1, N, sim.num, tn.diff) {
  ytn <- array(NA, dim = c(N,p.dim,sim.num))
  atn <- array(NA, dim = c(N,m.dim,sim.num))
  
  # simulate ytn
  for (i in 1:sim.num) {
    # simulate epsilon_tn
    epsilon.tn <- mvrnorm(n = N, mu = rep(0, p.dim), Sigma = H.matrix)
    
    # simulate eta_tn
    eta.tn <- matrix(NA, nrow = N, ncol = m.dim)
    for (tn in 1:(N-1)) {
      
      eta.tn[tn,] <- mvrnorm(n = 1, mu = rep(0, m.dim), Sigma = Q.matrix[[as.character(tn.diff[tn])]])
    }
    
    atn[1,,i] <- mvrnorm(n = 1, mu = rep(0, m.dim), Sigma = P1) # at0
    for (tn in 1:N) {
      ytn[tn,,i] <- mu + Z.matrix %*% atn[tn,,i] + epsilon.tn[tn,]
      if (tn < N) {
        atn[(tn + 1),,i] <- C.matrix[[as.character(tn.diff[tn])]] %*% atn[tn,,i] + eta.tn[tn,]
      }
    }
  }
  result <- list("ytn" = ytn, "atn" = atn)
  
  return(result)
}

### own neglogL
# neglogL.deter.fix <- function(para, deter, p.dim, m.dim, N, ytn, tn.diff) {
#   theta.diag <- exp(para[1:m.dim])
#   para <- para[-c(1:m.dim)]
#   Z.mat <- matrix(para[1:(p.dim*m.dim)], nrow = p.dim, byrow = TRUE)
#   para <- para[-c(1:(p.dim*m.dim))]
#   H.diag <- exp(para[1:p.dim])
#   para <- para[-c(1:p.dim)]
#   mu <- para
#   
#   C.matrix <- matrix(0, nrow = m.dim, ncol = m.dim)
#   theta.diag <- rev(cumsum(rev(theta.diag)))
#   if (deter > 0) {
#     theta.diag[1] <- theta.diag[1] + sqrt(deter)
#   }
#   theta.off.diag <- sqrt(((theta.diag[1]-theta.diag[2])^2 - deter) / 4)
#   C.matrix[upper.tri(C.matrix)] <- theta.off.diag
#   C.matrix[lower.tri(C.matrix)] <- -t(C.matrix)[lower.tri(-t(C.matrix))]
#   C.matrix <- C.matrix + diag(theta.diag)
#   C.matrix <- expm(-C.matrix)
#   H.matrix <- diag(H.diag) 
#   Q.matrix <- diag(1, m.dim)
#   
#   Ptn <- array(NA, dim = c(m.dim,m.dim,N))
#   Ptn[,,1] <- Q.matrix
#   atn <- matrix(0, nrow = N, ncol = m.dim)
#   
#   Ftn <- array(NA, dim = c(p.dim,p.dim,N))
#   Ftn[,,1] <- Z.mat %*% Ptn[,,1] %*% t(Z.mat) + H.matrix
# 
#   Ktn <- array(NA, dim = c(m.dim,p.dim,(N-1)))
#   Ktn[,,1] <- C.matrix %*% Ptn[,,1] %*% t(Z.mat) %*% solve(Ftn[,,1])
#   
#   vtn <- matrix(NA, nrow = N, ncol = p.dim)
#   vtn[1,] <- ytn[1,] - mu - Z.mat %*% atn[1,]
#   
#   # !!!!!! -(N + 1)*p/2 * log(2 * pi)
#   if (p.dim > 1) {
#     neg.logL <- (log(det(Ftn[,,1])) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
#   }
#   else {
#     neg.logL <- (log(Ftn[,,1]) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
#   }
#   
#   for (i in 1:(N-1)) {
#     atn[(i+1),] <- C.matrix %*% atn[i,] + Ktn[,,i] %*% vtn[i,]
#     Ptn[,,(i+1)] <- C.matrix %*% Ptn[,,i] %*% t(C.matrix - Ktn[,,i] %*% Z.mat) + Q.matrix
#     Ftn[,,(i+1)] <- Z.mat %*% Ptn[,,(i+1)] %*% t(Z.mat) + H.matrix
#     if (i < (N-1)) {
#       Ktn[,,(i+1)] <- C.matrix %*% Ptn[,,(i+1)] %*% t(Z.mat) %*% solve(Ftn[,,(i+1)])
#     }
#     vtn[(i+1),] <- ytn[(i+1),] - mu - Z.mat %*% atn[(i+1),]
#     if (p.dim > 1) {
#       neg.logL <- neg.logL + (log(det(Ftn[,,(i+1)])) + t(vtn[(i+1),]) %*% solve(Ftn[,,(i+1)]) %*% vtn[(i+1),]) / 2
#     }
#     else {
#       neg.logL <- neg.logL + (log(Ftn[,,(i+1)]) + t(vtn[(i+1),]) %*% solve(Ftn[,,(i+1)]) %*% vtn[(i+1),]) / 2
#     }
#     
#   }
#   
#   return(neg.logL)
# }


### own neglogL
# neglogL.Zfix <- function(para, Z.mat, p.dim, m.dim, N, ytn, tn.diff) {
#   theta.off.diag <- para[1:(m.dim*(m.dim-1)/2)]
#   para <- para[-c(1:(m.dim*(m.dim-1)/2))]
#   theta.diag <- exp(para[1:m.dim])
#   para <- para[-c(1:m.dim)]
#   H.diag <- exp(para[1:p.dim])
#   para <- para[-c(1:p.dim)]
#   mu <- para
#   
#   unique.t <- sort(unique(tn.diff))
#   theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
#   theta.mat[upper.tri(theta.mat)] <- theta.off.diag
#   theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
#   theta.diag <- rev(cumsum(rev(theta.diag)))
#   theta.mat <- theta.mat + diag(theta.diag)
#   
#   Z.mat <- matrix(Z.mat, nrow = p.dim, byrow = TRUE)
#   C.matrix.dic <- hash()
#   for (t in unique.t) {
#     C.matrix <- expm(-theta.mat*t)
#     C.matrix.dic[[as.character(t)]] <- C.matrix
#   }
#   
#   H.matrix <- diag(H.diag) 
#   Q.matrix.dic <- hash()
#   Q.matrix.dic[[as.character(1)]] <- diag(1,m.dim)
#   if (max(unique.t) > 1) {
#     for (t in 2:max(unique.t)) {
#       Q.matrix.dic[[as.character(t)]] <- Q.matrix.dic[[as.character(t - 1)]] + expm(-theta.mat*(t-1)) %*% t(expm(-theta.mat*(t-1)))
#     }
#   }
#   
#   Ptn <- array(NA, dim = c(m.dim,m.dim,N))
#   Ptn[,,1] <- diag(1,m.dim)
#   atn <- matrix(0, nrow = N, ncol = m.dim)
#   
#   Ftn <- array(NA, dim = c(p.dim,p.dim,N))
#   Ftn[,,1] <- Z.mat %*% Ptn[,,1] %*% t(Z.mat) + H.matrix
#   
#   Ktn <- array(NA, dim = c(m.dim,p.dim,(N-1)))
#   Ktn[,,1] <- C.matrix.dic[[as.character(1)]] %*% Ptn[,,1] %*% t(Z.mat) %*% solve(Ftn[,,1])
#   
#   vtn <- matrix(NA, nrow = N, ncol = p.dim)
#   vtn[1,] <- ytn[1,] - mu - Z.mat %*% atn[1,]
#   
#   # !!!!!! -(N + 1)*p/2 * log(2 * pi)
#   if (p.dim > 1) {
#     neg.logL <- (log(det(Ftn[,,1])) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
#   }
#   else {
#     neg.logL <- (log(Ftn[,,1]) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
#   }
#   
#   for (i in 1:(N-1)) {
#     atn[(i+1),] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% atn[i,] + Ktn[,,i] %*% vtn[i,]
#     Ptn[,,(i+1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% Ptn[,,i] %*% t(C.matrix.dic[[as.character(tn.diff[i])]] - Ktn[,,i] %*% Z.mat) + Q.matrix.dic[[as.character(tn.diff[i])]]
#     Ftn[,,(i+1)] <- Z.mat %*% Ptn[,,(i+1)] %*% t(Z.mat) + H.matrix
#     if (i < (N-1)) {
#       Ktn[,,(i+1)] <- C.matrix.dic[[as.character(tn.diff[i+1])]] %*% Ptn[,,(i+1)] %*% t(Z.mat) %*% solve(Ftn[,,(i+1)])
#     }
#     vtn[(i+1),] <- ytn[(i+1),] - mu - Z.mat %*% atn[(i+1),]
#     if (p.dim > 1) {
#       neg.logL <- neg.logL + (log(det(Ftn[,,(i+1)])) + t(vtn[(i+1),]) %*% solve(Ftn[,,(i+1)]) %*% vtn[(i+1),]) / 2
#     }
#     else {
#       neg.logL <- neg.logL + (log(Ftn[,,(i+1)]) + t(vtn[(i+1),]) %*% solve(Ftn[,,(i+1)]) %*% vtn[(i+1),]) / 2
#     }
#     
#   }
#   
#   return(neg.logL)
# }


# estimation.OU.single.Zfix <- function(theta.off.diag.init, theta.diag.init, ytn.sim, N, p.dim, m.dim, tn.diff, Z.mat) {
#   temp.sim <- as.matrix(ytn.sim)
#   
#   H.par.init <- diag(cov(temp.sim)) / 3
#   mu.init <- colMeans(temp.sim)
#   init.par <- c(theta.off.diag.init, log(theta.diag.init), H.par.init, mu.init)
#   
#   result.sim <- optim(init.par, Z.mat = Z.mat, p.dim = p.dim, m.dim = m.dim, N = N, ytn = temp.sim, tn.diff = tn.diff, neglogL.Zfix, method = "L-BFGS-B", hessian = TRUE)
#   
#   result.mat <- para.to.model.Zfix(result.sim$par, m.dim, p.dim)
#   
#   theta.matrix.dist <- result.mat$theta
#   C.matrix.dist <- result.mat$C.matrix
#   H.matrix.dist <- result.mat$H.mat
#   mu.dist <- result.mat$mu
#   Hessian.record <- result.sim$hessian
#   L.dist <- result.sim$value
#   result.para.dist <- result.sim$par
#   
#   result <- list("theta.matrix.dist" = theta.matrix.dist, "C.matrix.dist" = C.matrix.dist, "H.matrix.dist" = H.matrix.dist, "mu.dist" = mu.dist, "Hessian.record" = Hessian.record, "L.dist" = L.dist, "result.para.dist" = result.para.dist)
#   
#   return(result)
# }
# 
# para.to.model.Zfix <- function(para, m.dim, p.dim){
#   theta.off.diag <- para[1:(m.dim*(m.dim-1)/2)]
#   para <- para[-c(1:(m.dim*(m.dim-1)/2))]
#   theta.diag <- exp(para[1:m.dim])
#   para <- para[-c(1:m.dim)]
#   H.diag <- exp(para[1:p.dim])
#   para <- para[-c(1:p.dim)]
#   mu <- para
#   
#   theta.mat <- matrix(0, nrow = m.dim, ncol = m.dim)
#   theta.mat[upper.tri(theta.mat)] <- theta.off.diag
#   theta.mat[lower.tri(theta.mat)] <- -t(theta.mat)[lower.tri(-t(theta.mat))]
#   theta.diag <- rev(cumsum(rev(theta.diag)))
#   theta.mat <- theta.mat + diag(theta.diag)
#   
#   C.matrix <- expm(-theta.mat)
#   
#   H.matrix <- diag(H.diag) 
#   
#   result <- list("C.matrix" = C.matrix, "H.mat" = H.matrix, "mu" = mu, "theta" = theta.mat)
#   return(result)
# }

### own neglogL
### own neglogL
neglogL.m1 <- function(para, p.dim, N, ytn, tn.diff) {
  theta <- exp(para[1])
  para <- para[-c(1)]
  Z.mat <- matrix(para[1:p.dim], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:p.dim)]
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  mu <- para
  
  unique.t <- sort(unique(tn.diff))
  
  C.matrix.dic <- hash()
  for (t in unique.t) {
    C.matrix.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  H.matrix <- diag(H.diag) 
  Q.matrix.dic <- hash()
  for (t in unique.t) {
    Q.matrix.dic[[as.character(t)]] <- (1-exp(-2*theta*t)) / (2*theta)
  }
  
  Ptn <- rep(NA, N)
  Ptn[1] <- 1/(2*theta)
  atn <- rep(0, N)
  
  Ftn <- array(NA, dim = c(p.dim,p.dim,N))
  Ftn[,,1] <- Z.mat %*% Ptn[1] %*% t(Z.mat) + H.matrix
  
  Ktn <- array(NA, dim = c(1,p.dim,(N-1)))
  Ktn[,,1] <- C.matrix.dic[[as.character(1)]] %*% Ptn[1] %*% t(Z.mat) %*% solve(Ftn[,,1])
  
  vtn <- matrix(NA, nrow = N, ncol = p.dim)
  vtn[1,] <- ytn[1,] - mu - Z.mat %*% atn[1]
  
  # !!!!!! -(N + 1)*p/2 * log(2 * pi)
  if (p.dim > 1) {
    neg.logL <- (log(det(Ftn[,,1])) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
  }
  else {
    neg.logL <- (log(Ftn[,,1]) + t(vtn[1,]) %*% solve(Ftn[,,1]) %*% vtn[1,]) / 2
  }
  
  for (i in 1:(N-1)) {
    atn[(i+1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% atn[i] + Ktn[,,i] %*% vtn[i,]
    Ptn[(i+1)] <- C.matrix.dic[[as.character(tn.diff[i])]] %*% Ptn[i] %*% t(C.matrix.dic[[as.character(tn.diff[i])]] - Ktn[,,i] %*% Z.mat) + Q.matrix.dic[[as.character(tn.diff[i])]]
    Ftn[,,(i+1)] <- Z.mat %*% Ptn[(i+1)] %*% t(Z.mat) + H.matrix
    if (i < (N-1)) {
      Ktn[,,(i+1)] <- C.matrix.dic[[as.character(tn.diff[i+1])]] %*% Ptn[(i+1)] %*% t(Z.mat) %*% solve(Ftn[,,(i+1)])
    }
    vtn[(i+1),] <- ytn[(i+1),] - mu - Z.mat %*% atn[(i+1)]
    if (p.dim > 1) {
      neg.logL <- neg.logL + (log(det(Ftn[,,(i+1)])) + t(vtn[(i+1),]) %*% solve(Ftn[,,(i+1)]) %*% vtn[(i+1),]) / 2
    }
    else {
      neg.logL <- neg.logL + (log(Ftn[,,(i+1)]) + t(vtn[(i+1),]) %*% solve(Ftn[,,(i+1)]) %*% vtn[(i+1),]) / 2
    }
    
  }
  
  return(neg.logL)
}

estimation.OU.single.m1 <- function(theta.init, ytn.sim, N, p.dim, tn.diff) {
  temp.sim <- as.matrix(ytn.sim)
  
  M <- cov(temp.sim[1:(N-1),], temp.sim[2:N,])
  M <- eigen(M + t(M))
  Z.init <- M$vectors[,1]
  Z.init <- c(t(Z.init))
  H.par.init <- diag(cov(temp.sim)) / 3
  mu.init <- colMeans(temp.sim)
  init.par <- c(log(theta.init), Z.init, H.par.init, mu.init)
  
  result.sim <- optim(init.par, p.dim = p.dim, N = N, ytn = temp.sim, tn.diff = tn.diff, neglogL.m1, method = "L-BFGS-B", hessian = TRUE)
  result.mat <- para.to.model.m1(result.sim$par, p.dim)
  
  theta.matrix.dist <- result.mat$theta
  C.matrix.dist <- result.mat$C.matrix
  Z.matrix.dist <- result.mat$Z.mat[1,] / abs(result.mat$Z.mat[1,]) * result.mat$Z.mat
  H.matrix.dist <- result.mat$H.mat
  mu.dist <- result.mat$mu
  Hessian.record <- result.sim$hessian
  L.dist <- result.sim$value
  result.para.dist <- result.sim$par
  
  AIC.value <- 2 * L.dist + 2 * length(result.para.dist)
  BIC.value <- 2 * L.dist + log(N) * length(result.para.dist)
  
  result <- list("theta.matrix.dist" = theta.matrix.dist, "C.matrix.dist" = C.matrix.dist, "Z.matrix.dist" = Z.matrix.dist, "H.matrix.dist" = H.matrix.dist, "mu.dist" = mu.dist, "Hessian.record" = Hessian.record, 
                 "L.dist" = L.dist, "result.para.dist" = result.para.dist, "AIC" = AIC.value, "BIC" = BIC.value)
  
  return(result)
}


para.to.model.m1 <- function(para, p.dim){
  theta <- exp(para[1])
  para <- para[-c(1)]
  Z.mat <- matrix(para[1:p.dim], nrow = p.dim, byrow = TRUE)
  para <- para[-c(1:p.dim)]
  H.diag <- exp(para[1:p.dim])
  para <- para[-c(1:p.dim)]
  mu <- para
  
  C.matrix <- exp(-theta)
  
  H.matrix <- diag(H.diag) 
  
  result <- list("C.matrix" = C.matrix, "Z.mat" = Z.mat, "H.mat" = H.matrix, "mu" = mu, "theta" = theta)
  return(result)
}


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
























