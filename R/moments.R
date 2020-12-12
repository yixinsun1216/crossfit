#' Double Machine Learning Moments
#'
#' Moment functions for partial linear model and partial linear poisson model.
#'
#' @param theta scalar or vector of the treatment effect parameter(s).
#' @param Y a vector of the outcome variable.
#' @param D a vector or data.frame of the treatment variable(s) of interest.
#' @param m a vector of E[D|X=x] values, where X is a vector of controls.
#' @param s a vector of E[Y|X=x] values, where X is a vector of controls.


# ============================================================================
# moment functions for linear model
# ============================================================================
# finite nuisance parameter approach
#' @export
psi_plr <- function(theta, Y, D, m, s){
  return(1/nrow(D) * t(D - m) %*% (Y - s - D  %*% theta))
}

#' @export
psi_plr_grad <- function(theta, Y, D, m, s){
  return(-1/nrow(D) * t(D- m) %*% D)
}

#' @export
psi_plr_op <- function(theta, Y, D, m, s) {
  theta <- matrix(theta, nrow = dim(D)[2])
  N <- nrow(D)
  op <- (D - m) * as.vector(Y - s - D %*% theta)
  return((1 / N) * t(op) %*% op)
}

# concentrating out approach -------------------------------------
#' @export
psi_plr_conc <- function(theta, Y, D, m, s){
  return(1/nrow(D) * t(D - m) %*% (Y - s - (D - m) %*% theta))
}

#' @export
psi_plr_grad_conc <- function(theta, Y, D, m, s){
  return(-1/nrow(D) * t(D - m) %*% (D - m))
}

#' @export
psi_plr_op_conc <- function(theta, Y, D, m, s) {
  theta <- matrix(theta, nrow = dim(D)[2])
  N <- nrow(D)
  op <- (D - m) * as.vector(Y - s - (D - m) %*% theta)
  return((1 / N) * t(op) %*% op)
}

# ============================================================================
# moment functions for poisson
# ============================================================================
#' @export
psi_plpr <- function(theta, Y, D, m, s){
  theta <- matrix(theta, nrow = dim(D)[2])

  N <- nrow(D)
  return((1/N) * t(D - m) %*% (Y - exp(D %*% theta + s)))
}

#' @export
psi_plpr_grad <- function(theta, Y, D, m, s){
  s_list <- split(s, seq(nrow(s)))
  d_list <- split(D, seq(nrow(D)))
  N <- nrow(D)

  row_mult <- map2_dfr(d_list, s_list, function(d, s1) as.matrix(d) %*% exp(d %*% theta + s1))
  return(-1/N * t(D - m) %*% t(row_mult))
}

#' @export
psi_plpr_op <- function(theta, Y, D, m, s){
  theta <- matrix(theta, nrow = dim(D)[2])
  Z <- D - m
  z_list <- split(Z, seq(nrow(Z)))
  op_list <- split((Y - exp(D %*% theta + s)), seq(nrow(Y)))

  N <- nrow(D)
  op <- as.matrix(map2_dfr(z_list, op_list, function(x, y) x*y))
  return((1 / N) *  op %*% t(op))
}

# Concentrating out approach for poisson -------------------------------------
#' @export
# partially linear poisson model
# recall this is:
# (Y - exp(D*theta) * s(X) * z(theta, X)) *
# (D - exp(theta) * m(X) * z(theta, X))
# where z(theta, X) = 1 / (exp(theta) * m(X) + 1 - m(X))
psi_plpr_with_grad <- function(theta, Y, D, m, s) {
  # first verify that D is univariate and 0/1 valued
  K <- ncol(D)
  nvals <- length(unique(D))
  minval <- min(D)
  maxval <- max(D)

  if(!(K == 1 & nvals == 2 & minval == 0 & maxval == 1)) {
    stop("treatment variable not univariate and binary")
  } else {
    z <- 1 / (exp(theta) * m + 1 - m)
    A <- (Y - exp(D * theta) * s * z)
    B <- (D - exp(theta) * m * z)

    # psi is A * B so gradient is dA * B + A * dB
    dz <- -1 * (exp(theta) * m) * z^2
    dA <- -1 * (D * exp(D * theta) * s * z + exp(D * theta) * s * dz)
    dB <- -1 * (exp(theta) * m * z + exp(theta) * m * dz)

    return(list(A * B, dA * B + A * dB))
  }
}

#' @export
psi_plpr_conc <- function(theta, Y, D, m, s) {
  return(mean(psi_plpr_with_grad(theta, Y, D, m, s)[[1]]))
}

#' @export
psi_plpr_grad_conc <- function(theta, Y, D, m, s) {
  return(mean(psi_plpr_with_grad(theta, Y, D, m, s)[[2]]))
}

#' @export
psi_plpr_op_conc <- function(theta, Y, D, m, s) {
  vals <- psi_plpr_with_grad(theta, Y, D, m, s)[[1]]
  return(mean(vals^2))
}

# psi: function that gives the value of the Neyman-Orthogonal moment at a
  # given value of theta
# psi_grad: function that returns the gradient of psi with respect to theta
# psi_plr_op: function that gives the variance estimator at a given
  # value of theta.
