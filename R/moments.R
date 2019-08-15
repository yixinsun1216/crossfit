#' Double Machine Learning Moments
#'
#' Moment functions for partial linear model and partial linear poisson model.
#'
#' @param theta scalar or vector of the treatment effect parameter(s).
#' @param Y a vector of the outcome variable.
#' @param D a vector or data.frame of the treatment variable(s) of interest.
#' @param gamma a vector of E[Y|X=x], where X is a vector of controls.
#' @param delta a vector or data.frame of E[D|X=x].
#'


# partial linear model
# recall this is (Y - gamma(X) - theta * (D - delta(X))) * (D - delta(X))
psi_plr <- function(theta, Y, D, gamma, delta) {
  theta <- matrix(theta, nrow = dim(D)[2])
  N <- nrow(Y)
  return((1 / N) * t(D - delta) %*% (Y - gamma - (D - delta) %*% theta))
}

psi_plr_grad <- function(theta, Y, D, gamma, delta) {
  N <- nrow(Y)
  return(-1 * (1 / N) * t(D - delta) %*% (D - delta))
}

psi_plr_op <- function(theta, Y, D, gamma, delta) {
  theta <- matrix(theta, nrow = dim(D)[2])
  N <- nrow(Y)
  op <- (D - delta) * as.vector(Y - gamma - (D - delta) %*% theta)
  return((1 / N) * t(op) %*% op)
}

# partially linear poisson model
# recall this is:
# (Y - exp(D*theta) * gamma(X) * z(theta, X)) *
# (D - exp(theta) * delta(X) * z(theta, X))
# where z(theta, X) = 1 / (exp(theta) * delta(X) + 1 - delta(X))
psi_plpr_with_grad <- function(theta, Y, D, gamma, delta) {
  # first verify that D is univariate and 0/1 valued
  K <- ncol(D)
  nvals <- length(unique(D))
  minval <- min(D)
  maxval <- max(D)

  if(!(K == 1 & nvals == 2 & minval == 0 & maxval == 1)) {
    stop("treatment variable not univariate and binary")
  } else {
    z <- 1 / (exp(theta) * delta + 1 - delta)
    A <- (Y - exp(D * theta) * gamma * z)
    B <- (D - exp(theta) * delta * z)

    # psi is A * B so gradient is dA * B + A * dB
    dz <- -1 * (exp(theta) * delta) * z^2
    dA <- -1 * (D * exp(D * theta) * gamma * z + exp(D * theta) * gamma * dz)
    dB <- -1 * (exp(theta) * delta * z + exp(theta) * delta * dz)

    return(list(A * B, dA * B + A * dB))
  }
}

psi_plpr <- function(theta, Y, D, gamma, delta) {
  return(mean(psi_plpr_with_grad(theta, Y, D, gamma, delta)[[1]]))
}

psi_plpr_grad <- function(theta, Y, D, gamma, delta) {
  return(mean(psi_plpr_with_grad(theta, Y, D, gamma, delta)[[2]]))
}

psi_plpr_op <- function(theta, Y, D, gamma, delta) {
  vals <- psi_plpr_with_grad(theta, Y, D, gamma, delta)[[1]]
  return(mean(vals^2))
}
