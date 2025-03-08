% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.OU.single.R
\name{estimation.OU.single}
\alias{estimation.OU.single}
\title{Estimate Parameters for a Single Run of the Factor OUSSM Model}
\usage{
estimation.OU.single(
  theta.off.diag.init,
  theta.diag.init,
  ytn.sim,
  N,
  p.dim,
  m.dim,
  tn.diff
)
}
\arguments{
\item{theta.off.diag.init}{A numeric vector representing the initial values
for the off-diagonal elements of the state transition matrix \code{Theta}.}

\item{theta.diag.init}{A numeric vector representing the initial values for
the diagonal elements of \code{Theta}.}

\item{ytn.sim}{A matrix of dimensions \code{(N, p.dim)} containing
the observed data for estimation.}

\item{N}{An integer specifying the number of time points.}

\item{p.dim}{An integer specifying the number of observed variables.}

\item{m.dim}{An integer specifying the number of latent state variables.}

\item{tn.diff}{A numeric vector of length \code{N-1} representing the time
differences between consecutive observations.}
}
\value{
A list containing:
\item{theta.matrix.dist}{The estimated \code{Theta} matrix.}
\item{C.matrix.dist}{The estimated state transition matrix \code{C}.}
\item{Z.matrix.dist}{The estimated observation matrix \code{Z}.}
\item{H.matrix.dist}{The estimated observation noise covariance matrix \code{H}.}
\item{mu.dist}{The estimated mean vector.}
\item{Hessian.record}{The Hessian matrix from optimization.}
\item{L.dist}{The log-likelihood value from optimization.}
\item{result.para.dist}{The estimated parameter vector.}
\item{AIC}{The Akaike Information Criterion (AIC) value for model selection.}
\item{BIC}{The Bayesian Information Criterion (BIC) value for model selection.}
}
\description{
This function estimates the parameters of the Factor Ornstein-Uhlenbeck
State-Space Model (OUSSM) for a single dataset using Maximum Likelihood Estimation (MLE).
}
\details{
The function follows these steps:
\itemize{
\item Computes initial estimates for \code{Theta}, \code{Z}, and \code{H} based on the covariance of the data.
\item Uses eigen decomposition to initialize the observation matrix \code{Z}.
\item Runs the \code{optim()} function with the \code{L-BFGS-B} method to estimate parameters.
\item Converts the estimated parameters into structured matrices using \code{para.to.model()}.
\item Ensures sign consistency using \code{convert.Z.theta()}.
\item Computes model selection criteria (AIC and BIC).
}
}
\examples{
# Example usage (assuming required parameters and data are predefined)
result <- estimation.OU.single(theta.off.diag.init, theta.diag.init, ytn.sim, 
                               N, p.dim, m.dim, tn.diff)
result$theta.matrix.dist  # Access the estimated Theta matrix
}
