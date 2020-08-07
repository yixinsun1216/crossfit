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
    map(~dml_step_original(f, d, model, dml_seed, ml, poly_degree, family, ...), ...) %>%
    get_medians(nrow(d), dml_call)
}

dml_original_step <- function(f, d, model, ml, poly_degree, family){
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
  instruments <-
    map2(folds$train, folds$test,
         function(x, y) dml_fold(x, y, xnames, ynames, dnames, model, ml, family))
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
    list(D, Z) %>%
      pmap(function(D, Z) {
        psi_grad(D, Z)
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

  return(theta)
}

# =============================================================================
# Fill in components of moment conditions within each fold
# =============================================================================
# step 1. In training sample, perform linear lasso of y on x to select x_hat.
# Fit linear regression of y on x_hat and call these coefficients beta_hat
# step 2. Using test sample, calculate s = x_hat*beta_hat
# step 3. Using training data, perform linear lasso of d on x to select x_tilde.
# fit linear regressoin of d on x_tilde to get theta_tilde
# step 4. Calculate z = d - x_tilde*theta_tilde

dml_fold <- function(fold_train, fold_test, xnames, ynames, dnames, model, ml, family){

  # step 1 --------------------------------------------------------
  # formula for y on x
  f1 <-  paste(c(xnames), collapse = " + ") %>%
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

  # fit linear regression of y on x_hat to find beta_hat
  f2 <-  paste(x_hat, collapse = " + ") %>%
    paste0(ynames, " ~ ", .) %>%
    as.formula
  beta_hat <- tidy(lm(f2, fold_train))

  # step 2: s = x_hat*beta_hat -----------------------------
  s <- as.matrix(cbind(1, fold_test[x_hat])) %*% beta_hat$estimate

  # step 3 + 4 find theta_tilde & calculate z = d - x_tilde*theta_tilde ------
  # loop over each d to find theta_tilde
  z_k <- map_dfc(dnames, estimate_z, xnames, fold_train, fold_test, family)

  return(list(s = s, Z = z_k))
}

# Implement steps 3 & 4 -> calculate instrument z ------------------------------
estimate_z <- function(dvar, xnames, fold_train, fold_test, family){
  # formula for running d on x
  f_select <- paste(xnames, collapse = " + ") %>%
    paste0(dvar, " ~ ", .) %>%
    as.formula()

  # linear lasso of d on x to select x_tilde
  x_tilde <-
    coef(cv.glmnet(f_select, fold_train, family = family)) %>%
    tidy() %>%
    filter(row %in% xnames) %>%
    pull(row)

  # fit linear regression of d on x_tilde to find theta_tilde
  if(length(x_tilde) == 0){
    f_reg <- as.formula(paste(eval(dvar), "1", sep = "~"))
  } else{
    f_reg <- paste(x_tilde, collapse = " + ") %>%
      paste0(dvar, " ~ ", .) %>%
      as.formula
  }

  theta_tilde <- tidy(lm(f_reg, data = fold_train))$estimate

  # return z = d - x_tilde*theta_tilde
  d_test <- as.matrix(fold_test[[dvar]])
  x_test <- as.matrix(cbind(1, fold_test[x_tilde]))
  z_k <- tibble(!!dvar := as.vector(d_test - x_test %*% theta_tilde))
  return(z_k)
}


# ============================================================================
# moment functions for linear
# ============================================================================
psi_plr <- function(theta, Y, D, Z, s){
  return(1/length(Y) * t(Z) %*% (Y - s - D %*% theta))
}

psi_plr_grad <- function(D, Z){
  return(1/length(D)*t(Z) %*% D)
}


# ============================================================================
# moment functions for poisson
# ============================================================================
psi_plpr <- function(theta, Y, D, Z, s){
  theta <- matrix(theta, nrow = dim(D)[2])

  N <- nrow(Y)
  return((1/N) * t(Z) %*% (Y - exp(D %*% theta + s)))
}

# SIMPLIFY THIS FUNCTION
psi_plpr_grad <- function(theta, D, Z, s){
  s_list <- split(s, seq(nrow(s)))
  d_list <- split(D, seq(nrow(D)))
  z_list <- split(Z, seq(nrow(Z)))

  pmap(list(s_list, d_list, z_list), function(s1, d, z, ...) {
    -as.matrix(d) %*% exp(d %*% theta + s1) %*% z}) %>%
    reduce(`+`)/nrow(Z)
}

psi_plpr_op <- function(theta, Y, D, Z, s){
  theta <- matrix(theta, nrow = dim(D)[2])

  z_list <- split(Z, seq(nrow(Z)))
  op_list <- split((Y - exp(D %*% theta + s)), seq(nrow(Y)))

  N <- nrow(Y)
  op <- as.matrix(map2_dfr(z_list, op_list, function(x, y) x*y))
  return((1 / N) *  op %*% t(op))
}
