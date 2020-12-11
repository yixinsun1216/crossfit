#' Double Machine Learning Estimates
#'
#' @param f an object of class formula representing the model to be fitted.
#' @param d a dataframe containing the variables in f.
#' @param model model type or list of user created moment functions.
#'   The following model types are implementable: "linear" for partial linear
#'   model, "poisson" for a partial linear poisson model". If the argument is
#'   a list, the list must have three functions in order to generate theta,
#'   the coefficient of interest.
#' \enumerate{
#'   \item psi: function that gives the value of the Neyman-Orthogonal moment at a
#'   given value of theta
#'   \item psi_grad: function that returns the gradient of psi with respect to
#'   theta
#'   \item psi_plr_op: function that gives the variance estimator at a given
#'   value of theta.
#'   The default is `model = "linear"`.
#' @param ml Machine learning method to be used for estimating nuisance parameters.
#'   Currently it takes in `lasso`, and `rf` for regression forest. Note `rf`
#'   is not available for `score = "finite"`. Default is `ml = "lasso"`.
#' @param n Number of times to repeat the sample splitting and take median of
#'   results over the n samples. Default is `n = 100`.
#' @param k Number of folds for cross-fitting
#' @param score Takes either value `finite` or `concentrate`. `finite` refers to
#'   using the finite nuisance parameter orthogonal score construction, and
#'   `concentrate` refers to using the concentrating out approach.
#'   Default is `score = "finite"`
#' @param workers Number of workers to use in running the n dml calculations in
#'   parallel. Default is `workers = 1`, in which case the process is sequential.
#' @param drop_na if `TRUE`, then any row with an `NA` value is dropped. Default
#'   is `false`
#' @param family if `ml = lasso`, this is passed onto `cv.glmnet` to describe
#'    the response variable type.
#' @param poly_degree degree of polynomial for the nuisance parameters,
#'    to be used when `ml = "lasso"`. Default is `poly_degree = 1`.
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
#' Effect of temperature and precipitation on corn yield in the presence of
#' time and locational effects
#' data(corn_yield)
#' library(magrittr)
#'
#' yield_dml_rf <-
#'   "logcornyield ~ lower + higher + prec_lo + prec_hi | year + fips" %>%
#'   as.formula() %>%
#'   dml(corn_yield, "linear", n = 5, ml = "rf")
#'
#' yield_dml_lasso <-
#'   "logcornyield ~ lower + higher + prec_lo + prec_hi | year + fips" %>%
#'   as.formula() %>%
#'   dml(corn_yield, "linear", n = 5,  ml = "lasso", poly_degree = 3, score = "finite")
#'
#' @references V. Chernozhukov, D. Chetverikov, M. Demirer, E. Duflo, C. Hansen,
#' W. Newey, and J. Robins. Double/debiased machine learning for treatment and
#' structural parameters.The Econometrics Journal, 21(1):C1–C68, 2018a.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble enframe
#' @importFrom purrr pmap reduce map pluck map_dbl map2_dfr map_dfc map_lgl
#' @importFrom stats median optim update formula model.matrix model.frame
#' @importFrom furrr future_map future_options
#' @importFrom future plan multiprocess
#' @importFrom dplyr rename_all select bind_cols filter if_else mutate arrange pull
#' @importFrom stringr str_replace_all regex str_detect
#' @importFrom Formula Formula
#' @importFrom modelr crossv_kfold
#' @importFrom broom tidy
#' @importFrom glmnet cv.glmnet glmnet
#'
#'
#' @export
# main user-facing routine
#' @export
dml <- function(f, d, model = "linear", ml = "lasso", n = 101, k = 5,
                score = "concentrate", workers = 1, drop_na = FALSE,
                family = NULL, poly_degree = 1, lambda = NULL, ...) {
  dml_call <- match.call()

  # to make the median well-defined, add 1 to n if the user requests an even
  # number of splits
  nn <- if_else(n %% 2 == 0, n + 1, n)

  f <-
    Formula(f) %>%
    update(. ~ 0 + . | 0 + .)

  # lasso cannot handle NAs, so first prepare data if user specifies to drop NAs
  if(drop_na){
    d <- drop_na(d)
  }

  # parallelise the process if number of workers is > 1. Otherwise, future_map
  # defaults to sequential
  if(workers > 1){
    plan(future::multisession, workers = workers)
  }

  # make the estimation dataset -----------------------------------------------
  # (a) expand out any non-linear formula for y and sanitize names
  ty <- get_lhs_col(f, d)
  ynames <- names(ty)

  # (b) expand out any non-linear formula for d and sanitize names
  td <- get_rhs_cols(f, d, 1)
  dnames <- names(td)

  # (c) expand out any non-linear formula for x and sanitize names
  # expand out to polynomial depending on the user-inputted poly_degree
  tx <-
    get_rhs_cols(f, d, 2) %>%
    setNames(str_replace_all(names(.), "\\.|:", "\\_"))
  if(poly_degree > 1){
    tx <-
      poly(as.matrix(tx), degree = poly_degree) %>%
      as_tibble() %>%
      setNames(paste0("c", str_replace_all(names(.), "\\.|:", "\\_")))
  }
  xnames <- names(tx)

  # (d) make a new dataset of transformed y, transformed d and (transformed) x
  newdata <- bind_cols(ty, td, tx)

  # Calculate lambda to be used in cv.glmnet operation -------------------
  if(ml == "lasso" & is.null(lambda)){
    if(score == "finite"){
      y_reg <- cv.glmnet(as.matrix(cbind(tx, td)), as.matrix(ty),
                         standardize = FALSE, ...)
      l1 <- seq(y_reg$lambda.min, y_reg$lambda.1se, length.out = 10)
    }

    if(score == "concentrate"){
      y_reg <- cv.glmnet(as.matrix(tx), as.matrix(ty), standardize = FALSE, ...)
      l1 <- seq(y_reg$lambda.min, y_reg$lambda.1se, length.out = 10)
    }

    d_reg <-
      map(td, c) %>%
      map(function(d) cv.glmnet(as.matrix(tx), d, standardize = FALSE, ...))

    l2 <-
      d_reg %>%
      map(function(x) seq(x$lambda.min, x$lambda.1se, length.out = 10))

  } else{
    l1 <- lambda[[1]]
    l2 <- lambda[-1]
  }

  # pass into main dml function that is run n times -------------------
  seq(1, nn) %>%
    future_map(function(.x, ...)
      dml_step(f, newdata, tx, ty, td, model, ml, poly_degree,
               family, score, k, l1, l2, ...), ...,
               .options = future_options(packages = c("splines"))) %>%
    get_medians(nrow(d), dml_call)
}

dml_step <- function(f, newdata, tx, ty, td, model, ml, poly_degree,
                     family, score, k, l1, l2,  ...){
  # assign proper score function depending on if the user specifies using the
  # finite nuisance parameter vs concentrating-out approach, and whether the
  # model is linear or poisson
  if(model == "poisson" & score == "finite"){
    psi <- psi_plpr
    psi_grad <- psi_plpr_grad
    psi_op <- psi_plpr_op
  }
  if(model == "poisson" & score == "concentrate"){
    psi <- psi_plpr_conc
    psi_grad <- psi_plpr_grad_conc
    psi_op <- psi_plpr_op_conc
  }
  if(model == "linear" & score == "finite"){
    psi <- psi_plr
    psi_grad <- psi_plr_grad
    psi_op <- psi_plr_op
  }
  if(model == "linear" & score == "concentrate"){
    psi <- psi_plr_conc
    psi_grad <- psi_plr_grad_conc
    psi_op <- psi_plr_op_conc
  }
  if(score == "finite" & ml == "rf"){
    stop("rf only available for concentrating out method")
  }

  ynames <- names(ty)
  dnames <- names(td)
  xnames <- names(tx)

  # generate cross validation folds of this. Default is 5 folds
  folds <- crossv_kfold(newdata, k)
  folds$train <- map(folds$train, as_tibble)
  folds$test <- map(folds$test, as_tibble)

  # Calculate instruments -----------------------------------------------------
  # For each fold, calculate the partial outcome, s = x*beta, and m = x*mu
  if(score == "finite"){
    instruments <-
      map2(folds$train, folds$test, function(x, y)
        dml_fold(x, y, xnames, ynames, dnames, model, family, ml, l1, l2, ...))
  }

  # For the concentrating out approach, calculate the partial outcome, s = E[Y|X],
  # and m = E[D|X]
  if(score == "concentrate"){
    instruments <-
      map2(folds$train, folds$test, function(x, y)
        dml_fold_concentrate(x, y, xnames, ynames, dnames, ml, l1, l2, ...))
  }

  s <- map(instruments, pluck(1))
  m <- map(instruments, pluck(2))

  # Optimize over sample moment to find theta ---------------------------------
  Y <- map(folds$test, ~ as.matrix(select(., !!ynames)))
  D <- map(folds$test, ~ as.matrix(select(., !!dnames)))

  obj <- function(theta) {
    list(Y, D, m, s) %>%
      pmap(function(Y, D, m, s) {
        psi(theta, Y, D, m, s)
      }) %>%
      reduce(`+`) / length(D)
  }

  grad <- function(theta) {
    list(Y, D, m, s) %>%
      pmap(function(Y, D, m, s) {
        psi_grad(theta, Y, D, m, s)
      }) %>%
      reduce(`+`) / length(D)
  }

  theta0 <- rep(0, dim(D[[1]])[2])
  names(theta0) <- colnames(D[[1]])

  theta <-
    optim(theta0,
          function(x) square(obj(x)),
          function(x) grad(x) %*% obj(x),
          method = "BFGS",
          control = list(maxit = 5000))

  # Calculate covariance matrix -----------------------------------------------
  J0 <- grad(theta$par)
  s2 <-
    list(Y, D, m, s) %>%
    pmap(function(Y, D, m, s) {
      psi_op(theta$par, Y, D, m, s)
    }) %>%
    reduce(`+`) / length(D)

  s2 <- tryCatch(solve(t(J0), tol = 1e-20) %*% s2 %*% solve(J0),
                 error = function(e) return(NULL))

  return(list(theta$par, s2))
}

# =============================================================================
# Fill in components of moment conditions within each fold
# finite nuisance parameter approach
# =============================================================================
# step 1. In training sample, run lasso of y on x to select x_hat.
# for lasso, fit regression of y on x_hat and call these coefficients beta_hat
# step 2. Using test sample, calculate s = x_hat*beta_hat.
# calculate weights to be used for next round of lasso
# step 3. Using training data, perform weighted linear lasso of d on x to select
# x_tilde. For lasso, fit linear regression of d on x_tilde to get mu
# step 4. Calculate m = x_tilde*mu

dml_fold <- function(fold_train, fold_test, xnames, ynames, dnames, model,
                     family, ml, l1, l2, ...){
  # if using glmnet, make sure family is properly assigned
  # if model is linear, then use binomial for binary variable, and gaussian
  # otherwise
  if(is.null(family)){
    if(model == "poisson"){
      family <- "poisson"
    }else if(model == "linear" & length(unique(fold_train[,ynames])) == 2){
      family <- "binomial"
    } else{
      family <- "gaussian"
    }
  }
  # step 1 + 2 Linear --------------------------------------------------------
  # 1. use lasso of y on x to select x variables for linear model
  # 2. Calculate s = x_hat*beta_hat
  dep <- as.matrix(fold_train[, c(dnames, xnames)])
  resp <- as.matrix(fold_train[, ynames])

  x_hat <-
    coef(cv.glmnet(dep, resp, family = family, standardize = FALSE,
                   lambda = l1, ...)) %>%
    tidy() %>%
    filter(row %in% xnames) %>%
    pull(row)

  # if no x's are chosen, run regression of y on 1
  if(length(x_hat) == 0){
    f2 <-  paste(dnames, collapse = " + ") %>%
      paste0(ynames, " ~ ", .) %>%
      as.formula
  } else{
    f2 <-  paste(c(dnames, x_hat), collapse = " + ") %>%
      paste0(ynames, " ~ ", .) %>%
      as.formula
  }

  # fit regression of y on selected x
  beta_hat <- glm(f2, family, fold_train)
  coefs <-
    tidy(beta_hat) %>%
    filter(term %in% c("(Intercept)", x_hat))

  # Find s = x_hat*beta_hat
  s <- cbind(1, as.matrix(fold_test[,x_hat])) %*% as.matrix(coefs$estimate)

  if(model == "linear"){
    # weights for OLS is just 1s
    weights <- rep(1, nrow(fold_train))
  }
  if(model == "poisson"){
    # calculate weights w = exp(x_hat*beta_hat + d*theta_hat)
    weights <- exp(predict(beta_hat, data = fold_train))
  }

  # step 3 + 4 find mu & calculate m = x_tilde*mu -------------------------
  # loop over each d to find mu
  m_k <- map2_dfc(dnames, l2, function(x, y)
    estimate_m(dvar = x, l2 = y, xnames, weights, fold_train, fold_test, ml, ...))

  return(list(s = s, m = as.matrix(m_k)))
}

# =======================================================================
# Using concentrating out approach
# =======================================================================
# step 1. In training sample, train an ml model of y on x. Use this model to
# predict, for the testing sample s = E[Y|X]
# step 2. Using training data, train an ml model of d on x. Use this model to
# predict, for the test sample, n = E[Y|X]
dml_fold_concentrate <- function(fold_train, fold_test, xnames, ynames, dnames,
                                 ml, l1, l2,  ...){
  if(ml == "lasso"){
    # step 1:
    # pluck out x and y variables
    dep <- as.matrix(fold_train[, xnames])
    resp <- as.matrix(fold_train[, ynames])

    x_hat <-
      coef(cv.glmnet(dep, resp, lambda = l1, standardize = FALSE, ...)) %>%
      tidy() %>%
      filter(row %in% xnames) %>%
      pull(row)

    # if no x's are chosen, run regression of y on 1
    if(length(x_hat) == 0){
      f2 <- paste(ynames, "1", sep = " ~ ") %>%
        as.formula
    } else{
      f2 <-  paste(x_hat, collapse = " + ") %>%
        paste0(ynames, " ~ ", .) %>%
        as.formula
    }

    # fit regression of y on selected x
    beta_hat <- lm(f2, fold_train)

    # Find s = E[Y|X]
    s <- predict(beta_hat, fold_test)
  }

  if(ml == "rf"){
    # train model of Y on X
    f1 <-  paste(xnames, collapse = " + ") %>%
      paste0(ynames, " ~ ", .) %>%
      as.formula
    x_hat <- regression_forest2(f1, fold_train, ...)

    # Find s = E[Y|X]
    s <- as.matrix(predict_rf2(x_hat, fold_test))
  }

  weights <- rep(1, nrow(fold_train))

  # step 2 calculate m = E[D|X] -------------------------------------
  # loop over each d and concatenate results
  if(is.null(l2)) l2 <- rep(NA, length(dnames))
  m_k <- map2_dfc(dnames, l2, function(x, y)
    estimate_m(dvar = x, l2 = y, xnames, weights, fold_train, fold_test, ml, ...))

  return(list(s = s, m = as.matrix(m_k)))
}

# =======================================================================
# Function to calculate m = x_tilde*mu
# =======================================================================
estimate_m <- function(dvar, xnames, w, fold_train, fold_test, ml, l2, ...){
  # linear lasso of d on x to select x_tilde
  if(ml == "lasso"){
    # pluck out x and d
    dep <- as.matrix(fold_train[, xnames])
    resp <- as.matrix(fold_train[, dvar])

    x_tilde <-
      coef(cv.glmnet(dep, resp, weights = w, standardize = FALSE,
                     lambda = l2, ...)) %>%
      tidy() %>%
      filter(row %in% xnames) %>%
      pull(row)

    # fit linear regression of d on x_tilde to find mu
    if(length(x_tilde) == 0){
      f_reg <- as.formula(paste(eval(dvar), "1", sep = "~"))
    } else{
      f_reg <- paste(x_tilde, collapse = " + ") %>%
        paste0(dvar, " ~ ", .) %>%
        as.formula
    }
    mu <- lm(f_reg, data = fold_train)
    m_k <- predict(mu, fold_test)
  }

  if(ml == "rf"){
    # Train model of D on X
    f_select <- paste(xnames, collapse = " + ") %>%
      paste0(dvar, " ~ ", .) %>%
      as.formula()
    mu <- regression_forest2(f_select, fold_train, ...)

    # find E[D|X]
    m_k <- predict_rf2(mu, fold_test)
  }

  # return m = x_tilde*mu
  m_k <- setNames(tibble(as.vector(m_k)), dvar)
  return(m_k)
}


# =============================================================================
# Final routine to return median estimates
# =============================================================================
# Takes a set of DML sample split estimates and returns the median point
# estimate and the "median" covariance matrix.
# As suggested in the DML paper, the "median" covariance matrix is selected
# using the matrix operator norm, which is the highest svd of a matrix
get_medians <- function(estimates, n, dml_call) {

  # drop estimates where s2 is NULL
  s2_all <- map(estimates, ~pluck(., 2))
  keep <- !map_lgl(s2_all, is.null)
  estimates <- estimates[keep]

  if(sum(!keep) > 0) message(paste(sum(!keep), "J_0 matrices could not be inverted"))

  # if there are an even number of estimates, drop a random 1
  if(length(estimates) %% 2 == 0){
    estimates <- estimates[-sample(1:length(estimates), 1)]
  }

  median_theta <-
    estimates %>%
    map(~ pluck(., 1)) %>%
    reduce(rbind)

  if(length(estimates) > 1){
    median_theta <-  apply(median_theta, 2, median)
  }

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
