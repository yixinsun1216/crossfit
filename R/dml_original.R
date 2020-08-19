# Stata Poisson lasso  https://www.stata.com/manuals/lassoxpopoisson.pdf
# Stata Linear lasso https://www.stata.com/manuals/lassoxporegress.pdf

#' @export
dml_original <- function(f, d, model, n = 101, nw = 4, dml_seed = NULL, ml,
                poly_degree = 3, drop_na = FALSE, family = "gaussian", ...) {
  dml_call <- match.call()
  dml_seed <- ifelse(is.null(dml_seed), FALSE, as.integer(dml_seed))

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

  seq(1, nn) %>%
    map(~dml_step_original(f, d, model, ml, poly_degree, family, ...), ...) %>%
    get_medians(nrow(d), dml_call)
}

#' @export
dml_step_original <- function(f, d, model, ml, poly_degree, family, ...){
  if(model == "poisson"){
    psi <<- psi_plpr
    psi_grad <<- psi_plpr_grad
    psi_op <<- psi_plpr_op
  }
  if(model == "linear"){
    psi <<- psi_plr
    psi_grad <<- psi_plr_grad
    psi_op <<- psi_plr_op
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
  folds <- crossv_kfold(newdata)
  folds$train <- map(folds$train, as_tibble)
  folds$test <- map(folds$test, as_tibble)


  # Calculate instruments -----------------------------------------------------
  # For each fold, calculate the partial outcome, s = x*beta, and the instrument,
  # z = d - x*theta
  tic("instruments")
  instruments <-
    map2(folds$train, folds$test,
         function(x, y) dml_fold(x, y, xnames, ynames, dnames, model, ml, family))
  toc()
  s <- map(instruments, pluck(1))
  Z <- map(instruments, pluck(2))

  # Optimize over sample moment to find theta ---------------------------------
  Y <- map(folds$test, ~ as.matrix(select(., !!ynames)))
  D <- map(folds$test, ~ as.matrix(select(., !!dnames)))

  obj <- function(theta) {
    list(Y, D, Z, s) %>%
      pmap(function(Y, D, Z, s) {
        psi(theta, Y, D, Z, s)
      }) %>%
      reduce(`+`) / length(D)
  }

  grad <- function(theta) {
    list(D, Z, s) %>%
      pmap(function(D, Z, s) {
        psi_grad(theta, D, Z, s)
      }) %>%
      reduce(`+`) / length(D)
  }

  theta0 <- rep(0, dim(D[[1]])[2])
  names(theta0) <- colnames(D[[1]])

  tic("optimization")
  theta <-
    optim(theta0,
          function(x) square(obj(x)),
          function(x) grad(x) %*% obj(x),
          method = "BFGS",
          control = list(maxit = 500000))
  toc()

  # Calculate covariance matrix -----------------------------------------------
  J0 <- grad(theta$par)
  s2 <-
    list(Y, D, Z, s) %>%
    pmap(function(Y, D,Z, s) {
      psi_op(theta$par, Y, D, Z, s)
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

dml_fold <- function(fold_train, fold_test, xnames, ynames, dnames, model, ml, family){
  # if using glmnet, make sure family is properly assigned
  # if model is not poisson, let family default to binomial/gaussian
  if(model == "poisson" & is.null(family)){
    family <- "poisson"
  }

  # step 1 + 2 Linear --------------------------------------------------------
  # 1. use lasso of y on x to select x variables for linear model
  # 2. Calculate s = x_hat*beta_hat
  if(model == "linear" & ml %in% c("lasso", "rlasso")){
    f1 <-  paste(xnames, collapse = " + ") %>%
      paste0(ynames, " ~ ", .) %>%
      as.formula

    if(ml == "lasso"){
      # use linear lasso of y on x to select x variables
      x_hat <-
        coef(cv.glmnet(f1, fold_train, family = family)) %>%
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

    # if no x's are chosen, run regression of y on 1
    if(length(x_hat) == 0){
      f2 <- as.formula(paste(ynames, "1", sep = "~"))
    } else{
      f2 <-  paste(x_hat, collapse = " + ") %>%
        paste0(ynames, " ~ ", .) %>%
        as.formula
    }

    # fit regression of y on selected x
    beta_hat <- lm(f2, fold_train)

    # Find s = x_hat*beta_hat
    s <- predict(beta_hat, newdata = fold_test)

    # weights for OLS is just 1s
    weights <- rep(1, nrow(fold_train))
  }

  # Steps 1 & 2 Poisson --------------------------------------------------------
  # use lasso of y on d and x to select x variables for poisson model
  if(model == "poisson" & ml == "lasso"){
    f1 <-  paste(c(dnames, xnames), collapse = " + ") %>%
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

    # if no x's are chosen, run regression of y on 1
    if(length(x_hat) == 0){
      f2 <- as.formula(paste(ynames, "1", sep = "~"))
    } else{
      f2 <-  paste(x_hat, collapse = " + ") %>%
        paste0(ynames, " ~ ", .) %>%
        as.formula
    }

    # fit regression of y on selected x
    beta_hat <- glm(f2, "poisson", fold_train)
    coefs <-
      tidy(beta_hat) %>%
      filter(term %in% c("(Intercept)", x_hat))

    # Find s = x_hat*beta_hat
    s <- cbind(1, as.matrix(fold_test[,x_hat])) %*% as.matrix(coefs$estimate)

    # calculate weights w = exp(x_hat*beta_hat + d*theta_hat)
    weights <- exp(predict(beta_hat, data = fold_test))
  }

  if(ml == "regression_forest"){
    f1 <-  paste(xnames, collapse = " + ") %>%
      paste0(ynames, " ~ ", .) %>%
      as.formula

    x_hat <- regression_forest2(f1, fold_train)
    s <- as.matrix(predict_rf2(x_hat, fold_test))
    weights <- NULL
  }

  # step 3 + 4 find mu & calculate z = d - x_tilde*mu ------
  # loop over each d to find mu
  z_k <- map_dfc(dnames, estimate_z, xnames, weights, fold_train, fold_test, ml)

  return(list(s = s, Z = as.matrix(z_k)))
}

# Implement steps 3 & 4 -> calculate instrument z ------------------------------
estimate_z <- function(dvar, xnames, w, fold_train, fold_test, ml){
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
    x_theta <- predict(mu, fold_test)
  }

  if(ml == "regression_forest"){
    mu <- regression_forest2(f_select, fold_train)
    x_theta <- predict_rf2(mu, fold_test)
  }

  # return z = d - x_tilde*mu
  d_test <- fold_test[,dvar]
  z_k <- setNames(tibble(as.vector(d_test - x_theta)), dvar)
  return(z_k)
}


# ============================================================================
# moment functions for linear
# ============================================================================
#' @export
psi_plr <- function(theta, Y, D, Z, s){
  return(1/nrow(D) * t(Z) %*% (Y - s - D  %*% theta))
}

#' @export
psi_plr_grad <- function(theta, D, Z, s){
  return(-1/nrow(D) * t(Z) %*% D)
}

#' @export
psi_plr_op <- function(theta, Y, D, Z, s) {
  theta <- matrix(theta, nrow = dim(D)[2])
  N <- nrow(D)
  op <- Z * as.vector(Y - s - D %*% theta)
  return((1 / N) * t(op) %*% op)
}


# ============================================================================
# moment functions for poisson
# ============================================================================
#' @export
psi_plpr <- function(theta, Y, D, Z, s){
  theta <- matrix(theta, nrow = dim(D)[2])

  N <- nrow(D)
  return((1/N) * t(Z) %*% (Y - exp(D %*% theta + s)))
}

#' @export
psi_plpr_grad <- function(theta, D, Z, s){
  s_list <- split(s, seq(nrow(s)))
  d_list <- split(D, seq(nrow(D)))
  N <- nrow(D)

  row_mult <- map2_dfr(d_list, s_list, function(d, s1) as.matrix(d) %*% exp(d %*% theta + s1))
  return(-1/N * t(Z) %*% t(row_mult))
}

#' @export
psi_plpr_op <- function(theta, Y, D, Z, s){
  theta <- matrix(theta, nrow = dim(D)[2])

  z_list <- split(Z, seq(nrow(Z)))
  op_list <- split((Y - exp(D %*% theta + s)), seq(nrow(Y)))

  N <- nrow(D)
  op <- as.matrix(map2_dfr(z_list, op_list, function(x, y) x*y))
  return((1 / N) *  op %*% t(op))
}


#
# # house this here for now as a dup from ml_functions.R - make sure to get rid of this once we combine the two methods
# predict_rf2 <- function(forest, newdata = NULL) {
#   f <- forest[["formula"]]
#
#   if(!is.null(newdata)) {
#     X <-
#       formula(f, rhs = 1, lhs = 0) %>%
#       update(~ 0 + .) %>%
#       model.matrix(newdata)
#     return(pluck(predict(forest, X), "predictions"))
#   } else {
#     return(pluck(predict(forest), "predictions"))
#   }
# }
