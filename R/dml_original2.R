# https://www.stata.com/manuals/lassoxpopoisson.pdf

#' @export
dml_original <- function(f, d, model, n = 101, nw = 4, dml_seed = NULL, ml,
                poly_degree = 1, drop_na = FALSE, family = "gaussian", ...) {
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

  # plan(future::multisession, .init = nw)
  # seq(1, nn) %>%
  #   future_map(function(.x, ...) dml_step_original(f, d, model, dml_seed, ml,
  #                                         poly_degree, ...), ...,
  #              .options = future_options(packages = c("splines"),
  #                                        seed = dml_seed)) %>%
  #   get_medians(nrow(d), dml_call)

  # seq(1, nn) %>%
  #   map(~dml_step_original(f, d, model, dml_seed, ml, poly_degree, family, ...), ...) %>%
  #   get_medians(nrow(d), dml_call)

  dml_step_original(f, d, model, dml_seed, ml, poly_degree, family, ...)
}



dml_step_original <- function(f, d, model, dml_seed = NULL,
                              ml, poly_degree = 1, family, ...){
  if(model == "poisson"){
    psi <<- psi_plpr
    psi_grad <<- psi_plpr_grad
    psi_op <<- psi_plpr_op
  }

  # step 0: make the estimation dataset
  # step 1. Calculate theta and beta using training data
  # step 2. Use training data and estimates from 1 to find mu
  # step 3. Optimize over psi using hold out data

  # step 0 ----------------------------
  # (a) expand out any non-linear formula for y and sanitize names
  ty <- get_lhs_col(f, d)
  ynames <- names(ty)

  # (b) expand out any non-linear formula for d and sanitize names
  td <- get_rhs_cols(f, d, 1)
  dnames <- names(td)

  # (c) expand out any non-linear formula for x and sanitize names
  # expand out to polynomial depending on the user-inputted poly_degree
  # NOTE: no real reason to have nonlinear formulae here but might as well
  # be robust to it
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
  folds$train <- map(folds$train, data.frame)
  folds$test <- map(folds$test, data.frame)

  # Step 1  ----------------------------
  # 2 things we want to do:
  # i) Y = f(D*theta + X*beta, noise): Run poisson lasso of Y on X and D (
    # lasso only on X's) to pick out important X's
  # ii) D = g(X*gamma, noise): run normal poisson to estimate theta_hat and
    # beta_hat with all of D and the chosen subset of X from (i)
  # estimate_weights() returns beta_hat, as well as weights to be
    # used in the next step
  coefs <- estimate_weights(xnames, dnames, ynames, folds, ml, model, family)

  s <- map(coefs, function(x) x[[1]])

  # Step 2 ----------------------------
  # Use estimated coefficients to find mu, and calculate Z = D - X*mu_hat
  Z <-
    pmap(list(coefs, folds$train, folds$test),
          function(x, y, z, ...) estimate_z(x[[2]], y, z, dnames, xnames, ml, model)) %>%
    map(as.matrix)

  # Step 3  ----------------------------
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
  theta <-
    optim(theta0,
          function(x) square(obj(x)),
          function(x) 2 * grad(x) %*% obj(x),
          method = "BFGS",
          control = list(maxit = 500))

  print(theta$value)

  # calculate covariance matrix
  J0 <- grad(theta$par)
  s2 <-
    list(Y, D, Z, s) %>%
    pmap(function(Y, D,Z, s) {
      psi_op(theta$par, Y, D, Z, s)
    }) %>%
    reduce(`+`) / length(D)

  #s2 <- solve(t(J0), tol = 1e-20) %*% s2 %*% solve(J0)

  return(list(theta$par, s2))
}


# estimate theta and beta of Y = D*theta + X*beta using ML ---------------------
# Step 1: Poisson lasso of Y on D and X to select controls, x_hat
# Step 2: Fit Poisson regression of y on D and X_hat
# Step 3: use the beta_hat and X's from estimation sample to calculate s = X*beta_hat
estimate_weights <- function(xnames, dnames, ynames, folds, ml, model, family, ...){
  f1 <-  paste(c(xnames, dnames), collapse = " + ") %>%
    paste0(ynames, " ~ ", .) %>%
    as.formula

  if(ml == "lasso"){
    x_hat <- map(folds$train, function(y) {
      # include all D in lasso
      include <- names(y) %in% dnames

      ml_coef <- coef(cv.glmnet(f1, y, family = family,  penalty.factor = include))
      data.frame(coef = pluck(dimnames(ml_coef), 1), value = matrix(ml_coef)) %>%
        mutate(coef = as.character(coef)) %>%
        filter(value != 0,
               coef %in% xnames) %>%
        pull(coef)
    })
  }


  # FIGURE OUT HOW TO FIX SO THAT D IS ALWAYS INCLUDED IN LASS0 ================
  if(ml == "rlasso"){
    x_hat <- map(folds$train, function(y) {
      # include all D in lasso
      include <- names(y) %in% dnames

      # find coefficients on lasso and return as dataframe
      coef(rlasso2(f1, y, penalty = include)) %>%
        data.frame(coef = names(.), value = .)%>%
        mutate(coef = as.character(coef)) %>%
        filter(value != 0,
               !str_detect(coef, regex("Intercept", ignore_case = TRUE)),
               !(coef %in% dnames)) %>%
        pull(coef)
    })
  }

  # run regression of y on d and x
  if(model == "poisson"){
    coef_hat <- map2(folds$train, x_hat, function(x, y){
       f2 <-  paste(c(dnames, y), collapse = " + ") %>%
         paste0(ynames, " ~ ", .) %>%
         as.formula
       glm(f2, data = x) %>%
         tidy() %>%
         select(term, estimate)
     })

    # calculate weights using the training sample
    weights <- map2(folds$train, coef_hat, function(x, y){
      vars <- y$term[-1]
      exp(cbind(1, as.matrix(x[vars])) %*% as.matrix(y["estimate"]))
    })

    beta_hat <- map(coef_hat, function(x)
      filter(x, term %in% xnames | str_detect(term, regex("intercept", ignore_case = TRUE))))
  }

  s <- map2(folds$test, beta_hat, function(x, y) {
    vars <- y$term[-1]
    cbind(1, as.matrix(x[vars])) %*% as.matrix(y$estimate)
  })

  output <- map2(s, weights, list)
  return(output)
}

# For each component of D -------------------------------------------------
# Step 1: perform linear lasso of d on x using weights from last step
# Step 2: fit a weighted OLS using of d on x, and return Z = D - X*mu
estimate_z <- function(w, train_fold, test_fold, dnames, xnames, ml, model){
  # Pull out all x covariates and d variables
  x_train <- as.matrix(train_fold[, xnames])
  d_train <- as.list(train_fold[, dnames])

  x_test <- as.matrix(test_fold[, xnames])
  d_test <- as.list(test_fold[, dnames])


  z_k <- tibble(.rows = nrow(test_fold))
  if(ml == "lasso"){
    for(i in 1:length(dnames)){
      # Step 1: perform linear lasso of d on x using weights
      est <- coef(cv.glmnet(x_train, d_train[[i]], weights = w))

      # extract covariates selected using the linear lasso
      covs <-
        data.frame(coef = pluck(dimnames(est), 1), value = matrix(est)) %>%
        mutate(coef = as.character(coef)) %>%
        filter(value != 0, coef %in% xnames) %>%
        pull(coef)

      # Step 2: fit a weighted OLS using of d on x
      if(length(covs) == 0){
        f_mu <- as.formula(paste(eval(dnames[i]), "1", sep = "~"))
      } else{
         f_mu <- paste(covs, collapse = " + ") %>%
        paste0(dnames[[i]], " ~ ", .) %>%
        as.formula
      }

      df <-
        as_tibble(x_train) %>%
        mutate(!!dnames[i] := d_train[[i]])

      mu <- as.matrix(tidy(lm(f_mu, data = df, weights = w))$estimate)

      # return Z = D - X*mu
      z_calc <- as.matrix(d_test[[i]]) - as.matrix(cbind(1, x_test[,covs])) %*% mu
      z_k <- mutate(z_k, !!dnames[i] := as.vector(z_calc))
    }
  }
  return(z_k)
}



# ============================================================================
# moment functions for poisson
# ============================================================================
psi_plpr <- function(theta, Y, D, Z, s){
  theta <- matrix(theta, nrow = dim(D)[2])

  N <- nrow(Y)
  return((1/N) * t(Z) %*% (Y - exp(D %*% theta + s)))
}

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
