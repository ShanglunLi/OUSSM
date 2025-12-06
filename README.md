# OUSSM: Factor Ornstein-Uhlenbeck State Space Models

The **OUSSM** package implements the framework described in **"Factor State Space Modelling of the Ornstein-Uhlenbeck Process with Measurement Error and its Application"** by Shanglun Li, Toby Kenney, and Hong Gu.

Unlike standard univariate models, this package implements a **Factor State Space Model**. It models high-dimensional observed time series (e.g., microbiome counts or spatial temperature grids) using a lower-dimensional latent OU process. It explicitly accounts for measurement noise and solves identifiability issues through parameter constraints.

## Key Features

* **Factor Modeling:** Modeled as $y_t = \mu + Z x_t + \epsilon_t$, where $y_t$ is $p$-dimensional and $x_t$ is $m$-dimensional ($m \le p$).
* **Latent Dynamics:** The latent factors $x_t$ follow a multivariate Ornstein-Uhlenbeck process with mean reversion matrix $\Theta$.
* **Identifiability:** Implements constraints on $\Theta$ (sum of antisymmetric and diagonal matrices) and $Z$ to ensure unique parameter estimation.
* **Block Diagonalization:** Includes tools to transform estimated matrices into interpretable block-diagonal forms (separating real and complex eigenvalues).

## Installation

You can install the development version of OUSSM from [GitHub](https://github.com/ShanglunLi/OUSSM) with:

``` r
# install.packages("devtools")
devtools::install_github("ShanglunLi/OUSSM")
```

## Quick Start Example

This example demonstrates how to simulate a 2-dimensional latent factor process driving a 4-dimensional observed time series, and then estimate the parameters using the package's core function.

``` r
library(OUSSM)

# ---------------------------------------------------------
# 1. Simulate Synthetic Data (Factor OUSSM)
# ---------------------------------------------------------
set.seed(123)
N <- 300          # Time points
p.dim <- 4        # Dimension of observed data
m.dim <- 2        # Dimension of latent factors
dt <- 1           # Time step

# True Latent Dynamics (AR(1) approx of OU)
# We use a simple decay matrix for simulation purposes
C_true <- matrix(c(0.8, 0.1, -0.1, 0.7), nrow=2) 
Z_true <- matrix(rnorm(p.dim * m.dim), nrow=p.dim)
H_true <- diag(rep(0.1, p.dim)) # Measurement noise
Q_true <- diag(c(0.5, 0.5))     # Process noise

# Generate Data
x <- matrix(0, nrow=N, ncol=m.dim)
y <- matrix(0, nrow=N, ncol=p.dim)

for(t in 2:N) {
  x[t,] <- C_true %*% x[t-1,] + rnorm(m.dim, 0, sqrt(diag(Q_true)))
  y[t,] <- Z_true %*% x[t,] + rnorm(p.dim, 0, sqrt(diag(H_true)))
}

tn.diff <- rep(dt, N-1) # Regular time intervals

# ---------------------------------------------------------
# 2. Estimate Parameters
# ---------------------------------------------------------
# The estimation function requires initial guesses for Theta components.
# Theta is parameterized as: Antisymmetric + Diagonal (decreasing).

# Initial guess for off-diagonal theta (m*(m-1)/2 elements)
theta.off.init <- rep(0.1, m.dim * (m.dim - 1) / 2)

# Initial guess for diagonal theta (m elements)
theta.diag.init <- c(1.0, 0.5) 

# Run Estimation
# Note: Z and H are initialized automatically based on covariance PCA
results <- estimation.OU.single(
  theta.off.diag.init = theta.off.init,
  theta.diag.init = theta.diag.init,
  ytn.sim = y,
  N = N,
  p.dim = p.dim,
  m.dim = m.dim,
  tn.diff = tn.diff
)

# ---------------------------------------------------------
# 3. View Results
# ---------------------------------------------------------
print("Estimated Theta Matrix:")
print(results$theta.matrix.dist)

print("Estimated Loading Matrix (Z):")
print(results$Z.matrix.dist)

# Check model fit criteria
print(paste("AIC:", results$AIC))
print(paste("BIC:", results$BIC))
```

## Post-Processing: Block Diagonalization

The package provides a helper to interpret the dynamics by block-diagonalizing the $\Theta$ matrix. This separates independent mean-reverting components from coupled oscillatory components.

```r
# Apply block diagonalization to the estimated results
diag_results <- block_diagonalize(
  Theta = results$theta.matrix.dist,
  Z = results$Z.matrix.dist,
  Sigma = diag(m.dim) # Assuming identity for process noise scaling
)

print("Block-Diagonalized Theta:")
print(diag_results$block_matrix)
```

## Citation

If you use this package, please cite the accompanying article:

> Li, S., Kenney, T., & Gu, H. (2025). Factor State Space Modelling of the Ornstein-Uhlenbeck Process with Measurement Error and its Application. 
