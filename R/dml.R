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
#'   Currently it takes in `lasso`, `rlasso` (from the `hdm` package), and `rf`
#'   for regression forest. Note `rlasso` is not available for `model = "poisson"`.
#'   Default is `ml = "lasso"`.
#' @param n Number of times to repeat the sample splitting and take median of
#'   results over the n samples. Default is `n = 100`.
#' @param nfolds Number of folds for cross-fitting
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
#' @param poly_degree degree of polynomial for the nuisance parameters.
#'   Default is `poly_degree = 1`.
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
#'   dml(corn_yield, "linear", n = 5,
#'   ml = "regression_forest")
#'
#' yield_dml_lasso <-
#'   "logcornyield ~ lower + higher + prec_lo + prec_hi | year + fips" %>%
#'   as.formula() %>%
#'   dml(corn_yield, "linear", n = 5,
#'   ml = "lasso", family = "gaussian")
#'
#' @references V. Chernozhukov, D. Chetverikov, M. Demirer, E. Duflo, C. Hansen,
#' W. Newey, and J. Robins. Double/debiased machine learning for treatment and
#' structural parameters.The Econometrics Journal, 21(1):C1–C68, 2018a.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble enframe
#' @importFrom purrr pmap reduce map pluck map_dbl map2_dfr map_dfc
#' @importFrom stats median optim update formula model.matrix model.frame
#' @importFrom furrr future_map future_options
#' @importFrom future plan multiprocess
#' @importFrom dplyr rename_all select bind_cols filter if_else mutate arrange pull
#' @importFrom stringr str_replace_all regex str_detect
#' @importFrom Formula Formula
#' @importFrom modelr crossv_kfold
#' @importFrom broom tidy
#' @importFrom glmnetUtils cv.glmnet
#' @importFrom hdm rlasso
#'
#'
#' @export
# main user-facing routine
#' @export
dml <- function(f, d, model = "linear", ml = "lasso", n = 101, nfolds = 5,
                score = "finite", workers = 1, drop_na = FALSE, family = NULL,
                poly_degree = 3,
                ...) {
  dml_call <- match.call()

  # to make the median well-defined, add 1 to n if the user requests an even
  # number of splits
  nn <- if_else(n %% 2 == 0, n + 1, n)

  f <-
    Formula(f) %>%
    update(. ~ 0 + . | 0 + .)

  # lasso cannot handle NAs, so first prepare data if user specifies to drop NAs
  if(drop_na){
    d <-
      get_all_vars(f, d) %>%
      filter(complete.cases(.)) %>%
      as_tibble
  }

  # parallelise the process if number of workers is > 1. Otherwise, future_map
  # defaults to sequential
  if(workers > 1){
    plan(future::multisession, workers = workers)
  }

  seq(1, nn) %>%
    future_map(function(.x, ...) dml_step(f, d, model, ml, poly_degree, family,
                                          score, nfolds, ...), ...,
               .options = future_options(packages = c("splines"))) %>%
    get_medians(nrow(d), dml_call)


}

#' @export
dml_step <- function(f, d, model, ml, poly_degree, family, score, nfolds, ...){
  if(model == "poisson" & score == "finite"){
    psi <<- psi_plpr
    psi_grad <<- psi_plpr_grad
    psi_op <<- psi_plpr_op
  }
  if(model == "poisson" & score == "concentrate"){
    psi <<- psi_plpr_conc
    psi_grad <<- psi_plpr_grad_conc
    psi_op <<- psi_plpr_op_conc
  }
  if(model == "linear" & score == "finite"){
    psi <<- psi_plr
    psi_grad <<- psi_plr_grad
    psi_op <<- psi_plr_op
  }
  if(model == "linear" & score == "concentrate"){
    psi <<- psi_plr_conc
    psi_grad <<- psi_plr_grad_conc
    psi_op <<- psi_plr_op_conc
  }
  if(is.list(model)){
    psi <- model[[1]]
    psi_grad <- model[[2]]
    psi_op <- model[[3]]
  }
  if(model == "poisson" & ml == "rlasso"){
    stop("rlasso currently only availalbe for linear model")
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
    as.matrix(get_rhs_cols(f, d, 2)) %>%
    poly(degree = poly_degree, raw = TRUE) %>%
    as_tibble() %>%
    setNames(paste0("c", str_replace_all(names(.), "\\.", "\\_")))
  xnames <- names(tx)

  # (d) make a new dataset of transformed y, transformed d and (transformed) x
  newdata <- bind_cols(ty, td, tx)

  # (e) finally, generate cross validation folds of this. Default is 5 folds
  folds <- crossv_kfold(newdata, nfolds)
  folds$train <- map(folds$train, as_tibble)
  folds$test <- map(folds$test, as_tibble)


  # Calculate instruments -----------------------------------------------------
  # For each fold, calculate the partial outcome, s = x*beta, and the instrument,
  # z = d - x*theta
  if(score == "finite"){
    instruments <-
      map2(folds$train, folds$test,
           function(x, y) dml_fold(x, y, xnames, ynames, dnames, model, ml, family))
  }
  if(score == "concentrate"){
    instruments <-
      map2(folds$train, folds$test, function(x, y)
        dml_fold_concentrate(x, y, xnames, ynames, dnames, model, ml, family))
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
          control = list(maxit = 500000))

  # Calculate covariance matrix -----------------------------------------------
  J0 <- grad(theta$par)
  s2 <-
    list(Y, D, m, s) %>%
    pmap(function(Y, D, m, s) {
      psi_op(theta$par, Y, D, m, s)
    }) %>%
    reduce(`+`) / length(D)

  s2 <- solve(t(J0), tol = 1e-20) %*% s2 %*% solve(J0)

  return(list(theta$par, s2))
}

# =============================================================================
# Fill in components of moment conditions within each fold
# =============================================================================
# step 1. In training sample, run ml of y on x to select x_hat.
# for lasso, fit regression of y on x_hat and call these coefficients beta_hat
# step 2. Using test sample, calculate s = x_hat*beta_hat.
# calculate weights to be used for next round of lasso
# step 3. Using training data, perform linear lasso/RF of d on x to select x_tilde.
# For lasso, fit linear regression of d on x_tilde to get mu
# step 4. Calculate z = d - x_tilde*mu

dml_fold <- function(fold_train, fold_test, xnames, ynames, dnames, model, ml,
                     family, score){
  # if using glmnet, make sure family is properly assigned
  # if model is linear, then use binomial for binary variable, and gaussian
  # otherwise
  if(is.null(family)){
    if(model == "poisson"){
      family <- "poisson"
    }
    if(model == "linear" & length(unique(fold_train[,ynames])) == 2){
      family <- "binomial"
    }
    else{
      family <- "gaussian"
    }
  }

  # step 1 + 2 Linear --------------------------------------------------------
  # 1. use lasso of y on x to select x variables for linear model
  # 2. Calculate s = x_hat*beta_hat
  f1 <-
    paste0(c(dnames, xnames), collapse = " + ") %>%
    paste0(ynames, " ~ ", .) %>%
    as.formula

  if(ml == "lasso"){
    include <- as.numeric(names(fold_train) %in% dnames)

    x_hat <-
      coef(cv.glmnet(f1, fold_train, family = family, penalty.factor = include)) %>%
      tidy() %>%
      filter(row %in% xnames) %>%
      pull(row)
  }

  if(ml == "rlasso"){
    x_hat <-
      coef(rlasso(f1, fold_train)) %>%
      tidy() %>%
      filter(names %in% xnames, x != 0) %>%
      pull(names)
  }

  if( ml %in% c("lasso", "rlasso")){
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
      weights <- exp(predict(beta_hat, data = fold_test, type = "response"))
    }
  }

  if(ml == "rf"){
    f1 <-  paste(xnames, collapse = " + ") %>%
      paste0(ynames, " ~ ", .) %>%
      as.formula

    x_hat <- regression_forest2(f1, fold_train)
    s <- as.matrix(predict_rf2(x_hat, fold_test))
    weights <- NULL
  }

  # step 3 + 4 find mu & calculate z = d - x_tilde*mu ------
  # loop over each d to find mu
  m_k <- map_dfc(dnames, estimate_m, xnames, weights, fold_train, fold_test, ml)

  return(list(s = s, m = as.matrix(m_k)))
}

# =======================================================================
# Using concentrating out approach
# =======================================================================
dml_fold_concentrate <- function(fold_train, fold_test, xnames, ynames, dnames,
                                 model, ml, family, score){
  # if using glmnet, make sure family is properly assigned
  # if model is not poisson, let family default to binomial/gaussian
  if(is.null(family)){
    if(model == "poisson"){
      family <- "poisson"
    }
    if(model == "linear" & length(unique(fold_train[,ynames])) == 2){
      family <- "binomial"
    }
    else{
      family <- "gaussian"
    }
  }

  f1 <-  paste(xnames, collapse = " + ") %>%
    paste0(ynames, " ~ ", .) %>%
    as.formula

  if(ml == "lasso"){
    include <- as.numeric(names(fold_train) %in% dnames)

    x_hat <-
      coef(cv.glmnet(f1, fold_train, family = family, penalty.factor = include)) %>%
      tidy() %>%
      filter(row %in% xnames) %>%
      pull(row)
  }

  if(ml == "rlasso"){
    x_hat <-
      coef(rlasso(f1, fold_train)) %>%
      tidy() %>%
      filter(names %in% xnames, x != 0) %>%
      pull(names)
  }

  if(ml %in% c("lasso", "rlasso")){
    # if no x's are chosen, run regression of y on 1
    if(length(x_hat) == 0){
      f2 <- paste0(ynames, "1", sep = " ~ ") %>%
        as.formula
    } else{
      f2 <-  paste("0", paste(x_hat, collapse = " + "), sep = "+") %>%
        paste0(ynames, " ~ ", .) %>%
        as.formula
    }

    # fit regression of y on selected x
    beta_hat <- glm(f2, family, fold_train)

    # Find s = x_hat*beta_hat
    s <- predict(beta_hat, fold_test, type = "response")

    if(model == "linear"){
      # weights for OLS is just 1s
      weights <- rep(1, nrow(fold_train))
    }
    if(model == "poisson"){
      # calculate weights w = exp(x_hat*beta_hat + d*theta_hat)
      weights <- exp(predict(beta_hat, data = fold_test, type = "response"))
    }
  }

  if(ml == "rf"){
    x_hat <- regression_forest2(f1, fold_train)
    s <- as.matrix(predict_rf2(x_hat, fold_test))
    weights <- NULL
  }

  # step 3 + 4 find mu & calculate z = d - x_tilde*mu ------
  # loop over each d to find mu
  m_k <- map_dfc(dnames, estimate_m, xnames, weights, fold_train, fold_test, ml)

  return(list(s = s, m = as.matrix(m_k)))
}

# =======================================================================
# Steps 3 & 4: Calculate m = x_tilde*mu
# =======================================================================
estimate_m <- function(dvar, xnames, w, fold_train, fold_test, ml){

  # formula for running d on x
  f_select <- paste(xnames, collapse = " + ") %>%
    paste0(dvar, " ~ ", .) %>%
    as.formula()

  # linear lasso of d on x to select x_tilde
  if(ml == "lasso"){
    x_tilde <-
      coef(cv.glmnet(f_select, fold_train, weights = w)) %>%
      tidy() %>%
      filter(row %in% xnames) %>%
      pull(row)
  }

  if(ml == "rlasso"){
    x_tilde <-
      coef(rlasso(f_select, fold_train)) %>%
      tidy() %>%
      filter(names %in% xnames, x != 0) %>%
      pull(names)
  }

  if(ml %in% c("lasso", "rlasso")){
    # fit linear regression of d on x_tilde to find mu
    if(length(x_tilde) == 0){
      f_reg <- as.formula(paste(eval(dvar), "1", sep = "~"))
    } else{
      f_reg <- paste(x_tilde, collapse = " + ") %>%
        paste0(dvar, " ~ ", .) %>%
        as.formula
    }
    mu <- lm(f_reg, data = fold_train)
    x_mu <- predict(mu, fold_test)
  }

  if(ml == "rf"){
    mu <- regression_forest2(f_select, fold_train)
    x_mu <- predict_rf2(mu, fold_test)
  }

  # return m = x_tilde*mu
  m_k <- setNames(tibble(as.vector(x_mu)), dvar)
  return(m_k)
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

# Stata Poisson lasso  https://www.stata.com/manuals/lassoxpopoisson.pdf
# Stata Linear lasso https://www.stata.com/manuals/lassoxporegress.pdf

