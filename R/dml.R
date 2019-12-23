#' Double Machine Learning Estimates
#'
#' @param f an object of class formula representing the model to be fitted.
#' @param d a dataframe containing the variables in f.
#' @param psi function that gives the value of the Neyman-Orthogonal moment at a
#'   given value of theta
#' @param psi_grad function that gives the gradient of psi with respect to theta
#' @param psi_plr_op function that gives the variance estimator at a given
#'   value of theta.
#' @param n
#' @param nw
#' @param dml_seed
#' @param ml
#'
#' @return
#' \code{dml} returns an object of class "dml" with the following components:
#' \describe{
#'    \item{coefficients}{a named vector of coefficients.}
#'    \item{vcov}{variance-covariance matrix of the main parameters.}
#'    \item{nobs}{number of observations used}
#'    \item{call}{original function call with given arguments}
#' }
#'
#' @examples
#' # Effect of temperature and precipitation on corn yield in the presence of
#' # time and locational effects
#' data(corn_yield)
#' library(magrittr)
#'
#' yield_dml <-
#'   "logcornyield ~ lower + higher + prec_lo + prec_hi | year + fips" %>%
#'   as.formula() %>%
#'   dml(corn_yield, psi_plr, psi_plr_grad, psi_plr_op, n = 5,
#'   ml = "regression_forest", dml_seed = 123)
#'
#' @references V. Chernozhukov, D. Chetverikov, M. Demirer, E. Duflo, C. Hansen,
#' W. Newey, and J. Robins. Double/debiased machine learning for treatment and
#' structural parameters.The Econometrics Journal, 21(1):C1–C68, 2018a.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble enframe
#' @importFrom purrr pmap reduce map pluck map_dbl
#' @importFrom stats median optim update formula model.matrix model.frame
#' @importFrom furrr future_map future_options
#' @importFrom future plan multiprocess
#' @importFrom dplyr rename_all select bind_cols filter if_else
#' @importFrom stringr str_replace_all regex
#' @importFrom Formula Formula
#' @importFrom modelr crossv_kfold
#'
#'
#' @export
# main user-facing routine
dml <- function(f, d, psi, psi_grad, psi_op, n = 101, nw = 4,
                dml_seed = NULL, ml,  poly_degree = 3, drop_na = FALSE, ...) {
  dml_call <- match.call()
  dml_seed <- ifelse(is.null(dml_seed), FALSE, as.integer(dml_seed))

  # to make the median well-defined, add 1 to n if the user requests an even
  # number of splits
  nn <- if_else(n %% 2 == 0, n + 1, n)

  # ensure fm has no intercept in the parametric part
  f <-
    Formula(f) %>%
    update(. ~ 0 + . | 0 + .)

  # lasso cannot handle NAs, so first prepare data if user specifies to drop NAs
  if(drop_na){
    d <-
      model.frame(f, d) %>%
      filter(complete.cases(.)) %>%
      as_tibble
  }

  plan(future::multisession, .init = nw)
  seq(1, nn) %>%
    future_map(function(.x, ...) dml_step(f, d, psi, psi_grad, psi_op, dml_seed, ml,
                         poly_degree, ...), ...,
               .options = future_options(packages = c("splines"),
                                         seed = dml_seed)) %>%
    get_medians(nrow(d), dml_call)

}

square <- function(x) 0.5 * t(x) %*% x


dml_estimate <- function(Y, D, gamma, delta, psi, psi_grad, psi_op,
                         bounds = NULL) {
  obj <- function(theta) {
    list(Y, D, gamma, delta) %>%
      pmap(function(Y, D, gamma, delta) {
        psi(theta, Y, D, gamma, delta)
      }) %>%
      reduce(`+`) / length(D)
  }

  grad <- function(theta) {
    list(Y, D, gamma, delta) %>%
      pmap(function(Y, D, gamma, delta) {
        psi_grad(theta, Y, D, gamma, delta)
      }) %>%
      reduce(`+`) / length(D)
  }

  theta0 <- rep(0, dim(D[[1]])[2])
  names(theta0) <- colnames(D[[1]])
  theta <-
    optim(theta0,
          function(x) square(obj(x)),
          function(x) 2 * grad(x) %*% obj(x),
          method = "BFGS",
          control = list(maxit = 500))

  J0 <-
    list(Y, D, gamma, delta) %>%
    pmap(function(Y, D, gamma, delta) {
      psi_grad(theta$par, Y, D, gamma, delta)
    }) %>%
    reduce(`+`) / length(D)

  s2 <-
    list(Y, D, gamma, delta) %>%
    pmap(function(Y, D, gamma, delta) {
      psi_op(theta$par, Y, D, gamma, delta)
    }) %>%
    reduce(`+`) / length(D)

  s2 <- solve(t(J0)) %*% s2 %*% solve(J0)

  # note that what we want for inference is actually s2/N, implemented in the
  # get_medians function below
  return(list(theta$par, s2))
}

get_rhs_cols <- function(f, d, part = 1) {
  options(na.action='na.pass')

  f %>%
    formula(lhs = 0, rhs = part) %>%
    model.matrix(d) %>%
    as_tibble %>%
    rename_all(~ str_replace_all(., regex("[(, =)]"), "_"))
}

get_lhs_col <- function(f, d) {
  f %>%
    formula(lhs = 1, rhs = 0) %>%
    model.frame(d, na.action = NULL) %>%
    as_tibble %>%
    rename_all(~ str_replace_all(., regex("[(, =)]"), "_"))
}

# estimate theta and s2 in a single sample split of the data
# formula should be y ~ d | x
dml_step <- function(f, d, psi, psi_grad, psi_op, dml_seed = NULL,
                     ml, poly_degree, ...) {
  # step 1: make the estimation dataset
  # (a) expand out any non-linear formula for y and sanitize names
  ty <- get_lhs_col(f, d)

  # (b) expand out any non-linear formula for d and sanitize names
  td <- get_rhs_cols(f, d, 1)

  # (c) expand out any non-linear formula for x and sanitize names
  # NOTE: no real reason to have nonlinear formulae here but might as well
  # be robust to it
  tx <- get_rhs_cols(f, d, 2)

  # (c) make a new dataset of transformed y, transformed d and (transformed) x
  newdata <- bind_cols(ty, td, tx)

  # (d) finally, generate cross validation folds of this
  folds <- crossv_kfold(newdata)

  # step 2: save formulae for gamma and delta, based on transformed y and d
  xvars <- paste("0", paste0(names(tx), collapse = " + "), sep = " + ")


  f_gamma <-
    paste(names(ty), xvars, sep = " ~ ") %>%
    as.formula

  fs_delta <-
    names(td) %>%
    map(~ paste(., xvars, sep = " ~ ")) %>%
    map(as.formula)

  # Next 2 steps done in the function run_ml()
  # step 3: train models for delta and gamma
  # step 4: estimate values of delta and gamma in the hold out sample

  output <- run_ml(f_gamma, fs_delta, folds, ml, dml_seed, poly_degree, ...)
  gamma <- output$gamma
  delta <- output$delta

  # step 5: put together the values of Y and D in the hold out sample, and pass
  # things off to the dml_estimate routine
  ynames <- names(ty)
  Y <- map(folds$test, ~ as.matrix(select(as.data.frame(.), !!ynames)))

  dnames <- names(td)
  D <- map(folds$test, ~ as.matrix(select(as.data.frame(.), !!dnames)))

  return(dml_estimate(Y, D, gamma, delta, psi, psi_grad, psi_op))
}


# final routine which takes a set of DML sample split estimates and returns the
# median point estimate and the "median" covariance matrix
# as suggested in the DML paper, the "median" covariance matrix is selected
# using the matrix operator norm, which is the highest svd of a matrix
get_medians <- function(estimates, n, dml_call) {
  median_theta <-
    estimates %>%
    map(~ pluck(., 1)) %>%
    reduce(rbind) %>%
    apply(2, median)

  names(median_theta) <- names(estimates[[1]][[1]])

  medsq <- function(x) (x - median_theta) %*% t(x - median_theta)
  s2s <-
    estimates %>%
    map(~ pluck(., 2) + medsq(pluck(., 1)))

  s2_norms <-
    s2s %>%
    map_dbl(~ norm(., type = "2")) %>%
    enframe %>%
    filter(value == median(value))

  median_s2 <- s2s[[s2_norms$name]]

  return(structure(list(coefficients = median_theta,
                        vcov = (1 / n) * median_s2,
                        nobs = n,
                        call = dml_call),
                   class = "dml"))
}


# Coefficients: tune = TRUE
#   Estimate Std. Error
# lower    2.927e-05      0.000
# higher  -2.049e-03      0.000
# prec_lo -3.912e-03      0.001
# prec_hi  2.662e-04      0.000
# 406 seconds, seed 123
# 20.2 seconds
