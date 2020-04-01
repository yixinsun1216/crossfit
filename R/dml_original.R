# https://www.stata.com/manuals/lassoxpopoisson.pdf

# Effect of temperature and precipitation on corn yield in the presence of
# time and locational effects
library(dplyr)
library(crossfit)
library(splines)
load("C:/Users/Yixin Sun/Dropbox (Personal)/EPIC/texas/generated_data/final_leases.Rda")

reg_data <-
  final_leases %>%
  filter(InSample) %>%
  mutate(Private = NParcels15 > 0 | Type == "RAL")

base_controls <-
  "Auction + bs(Acres, df = 7) + Term + RoyaltyRate"

bonus_formula <-
  paste("BonusPerAcre", base_controls, sep = " ~ ") %>%
  paste("CentLat + CentLong + EffDate", sep = " | ") %>%
  as.formula

check <- dml_original(bonus_formula, reg_data, "poisson", n = 4, nw = 1, dml_seed = NULL, ml = "lasso", poly_degree = 1, drop_na = FALSE)


#' @export
dml_original <- function(f, d, model, n = 101, nw = 4, dml_seed = NULL, ml,
                poly_degree = 1, drop_na = FALSE, ...) {
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

  seq(1, nn) %>%
    map(~dml_step_original(f, d, model, dml_seed, ml, poly_degree, ...), ...) %>%
    get_medians(nrow(d), dml_call)
}



dml_step_original <- function(f, d, model, dml_seed = NULL,
                              ml, poly_degree = 1, ...){
  if(model == "poisson"){
    psi <<- psi_plpr
    psi_grad <<- psi_plpr_grad
    psi_op <<- psi_plpr_op
  }

  # step 1: make the estimation dataset
  # (a) expand out any non-linear formula for y and sanitize names
  ty <- get_lhs_col(f, d)
  ynames <- names(ty)

  # (b) expand out any non-linear formula for d and sanitize names
  td <- get_rhs_cols(f, d, 1)
  dnames <- names(td)

  # (c) expand out any non-linear formula for x and sanitize names
  # NOTE: no real reason to have nonlinear formulae here but might as well
  # be robust to it
  tx <- get_rhs_cols(f, d, 2)
  xnames <- names(tx)

  # (c) make a new dataset of transformed y, transformed d and (transformed) x
  newdata <- bind_cols(ty, td, tx)

  # (d) finally, generate cross validation folds of this
  folds <- crossv_kfold(newdata)
  folds$train <- map(folds$train, data.frame)

  # 1. Calculate theta and beta using training data
  # 2. Use training data and estimates from 1 to find mu
  # 3. Optimize over psi using hold out data
  f_ml <- paste(c(xnames, dnames), collapse = " + ") %>%
    paste0(ynames, " ~ ", .) %>%
    as.formula

  # Step 1:
  ml_coefs <- estimate_ml(f_ml, folds, ml, model, ynames, ...)

  # Step 2: Use estimated coefficients to find
  # mu = E[X'exp(D*theta + X*beta)X]^{-1}E[X'exp(D*theta - X*beta)D]
  mu_hat <-
    map2(ml_coefs, folds$train,
          function(x, y) estimate_mu(x, y, dnames, xnames))

  beta_hat <-
    map(ml_coefs, function(x) {
      filter(x, coef %in% xnames) %>%
        arrange(order(match(coef, xnames))) %>%
        pull(value)})

  # Step 3:
  Y <- map(folds$test, ~ as.matrix(select(as.data.frame(.), !!ynames)))
  D <- map(folds$test, ~ as.matrix(select(as.data.frame(.), !!dnames)))
  X <- map(folds$test, ~ as.matrix(select(as.data.frame(.), !!xnames)))

  obj <- function(theta) {
    list(Y, D, X, mu_hat, beta_hat) %>%
      pmap(function(Y, D, X, mu_hat, beta_hat) {
        psi(theta, Y, D, X, mu_hat, beta_hat)
      }) %>%
      reduce(`+`) / length(D)
  }

  grad <- function(theta) {
    list(D, X, mu_hat, beta_hat) %>%
      pmap(function(D, X, mu_hat, beta_hat) {
        psi_grad(theta, D, X, mu_hat, beta_hat)
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

  J0 <- grad(theta$par)

  s2 <-
    list(Y, D, X, mu_hat, beta_hat) %>%
    pmap(function(Y, D, X, mu_hat, beta_hat) {
      psi_op(theta$par, Y, D, X, mu_hat, beta_hat)
    }) %>%
    reduce(`+`) / length(D)

  s2 <- solve(t(J0), tol = 1e-20) %*% s2 %*% solve(J0)

  return(list(theta$par, s2))
}


# ============================================================================
# Use Original DML approach (not the concentrating out approach)
# ============================================================================
# estimate theta and beta of Y = D*theta + X*beta using ML
estimate_ml <- function(f_ml, folds, ml, model, ynames, ...){
  one_ml <- function(y, train_ml, ...){
    ml_coef <- coef(get(train_ml)(f_ml, y, ...))
    data.frame(coef = pluck(dimnames(ml_coef), 1), value = matrix(ml_coef)) %>%
      mutate(coef = as.character(coef)) %>%
      mutate(coef = if_else(str_detect(coef, regex("Intercept", ignore.case = TRUE)),
                            "Intercept", coef))
  }

  if(model == "poisson"){ f_ml <- update(f_ml, log(.) ~ .)}

  if(ml == "lasso"){
    ml_coefs <- map(folds$train, function(y) {
      ml_coef <- coef(cv.glmnet(f_ml, y, family = "gaussian", ...))
      data.frame(coef = pluck(dimnames(ml_coef), 1), value = matrix(ml_coef)) %>%
        mutate(coef = as.character(coef)) %>%
        mutate(coef = if_else(str_detect(coef, regex("Intercept", ignore.case = TRUE)),
                              "Intercept", coef))
    })
  }

  if(ml == "rlasso"){
    ml_coefs <- map(folds$train, function(y) {
      coef(rlasso2(f_ml, y, ...)) %>%
        data.frame(coef = names(.), value = .)%>%
        mutate(coef = as.character(coef)) %>%
        mutate(coef = if_else(str_detect(coef, regex("Intercept", ignore.case = TRUE)),
                              "Intercept", coef))
    })
  }

  return(ml_coefs)
}

# calculate mu from the estimated theta and beta coefficients, where mu is
# mu = E[X'exp(D*theta + X*beta)X]^{-1}E[X'exp(D*theta - X*beta)D]
estimate_mu <- function(coefs, train_fold, d_names, x_names){
  # shape D vars and coefs
  d_vars <- train_fold[, d_names]
  d_coefs <-
    filter(coefs, coef %in% d_names) %>%
    arrange(order(match(coef, d_names)))

  # Shape X vars and coefs
  x_vars <- train_fold[, x_names]
  x_coefs <-
    filter(coefs, coef %in% x_names) %>%
    arrange(order(match(coef, x_names)))

  intercept <- 0
  if(sum(str_detect(coefs$coef, "Intercept")) > 0){
    intercept <- coefs$value[coefs$coef == "Intercept"]
  }

  # calculate J_ThetaBeta and J_BetaBeta
  # for each row in the d_vars and x_vars data, calculate
    # 1. x_vars*exp(d_vars*theta + x_vars*beta)d_vars
    # 2. x_vars*exp(d_vars*theta + x_vars*beta)x_vars
  # calculate mu by dividing the sample averages of 1 by 2
  x_list <- split(x_vars, seq(nrow(x_vars)))
  d_list <- split(d_vars, seq(nrow(d_vars)))
  j_tb <-
    map2(x_list, d_list, function(x, y)
          j_calc(x, y, d_coefs$value, x_coefs$value, intercept, "tb")) %>%
    reduce(`+`)/length(x_list)

  j_bb <-
    map2(x_list, d_list, function(x, y)
      j_calc(x, y, d_coefs$value, x_coefs$value, intercept, "bb")) %>%
    reduce(`+`)/length(x_list)

  return(j_tb %*% solve(j_bb))
}

j_calc <- function(X, D, theta, beta, int, type){
  V <- D
  if(type == "bb") V <- X
  t(as.matrix(V*exp(sum(D*theta) + sum(X*beta) + int))) %*% as.matrix(X)
}


# ============================================================================
# moment functions for poisson
# ============================================================================
# psi_plpr <- function(theta, Y, D, X, mu, beta){
#   theta <- matrix(theta, nrow = dim(D)[2])
#   beta <- matrix(beta, nrow = dim(X)[2])
#   N <- nrow(Y)
#   return((1/N) * t(D - X %*% t(mu)) %*% (Y - exp(D %*% theta + X %*% beta)))
# }

psi_plpr <- function(theta, Y, D, X, mu, beta){
  y_list <- split(Y, seq(nrow(Y)))
  x_list <- split(X, seq(nrow(X)))
  d_list <- split(D, seq(nrow(D)))

  list(y_list, x_list, d_list) %>%
    pmap(function(y, x, d) {
      t(d - x %*% t(mu))%*% (y - exp(d %*% theta + x %*% beta))}) %>%
    reduce(`+`)/nrow(Y)
}


psi_plpr_grad <- function(theta, D, X, mu, beta){
  x_list <- split(X, seq(nrow(X)))
  d_list <- split(D, seq(nrow(D)))

  map2(x_list, d_list, function(x, d) {
      (-t(d - x %*% t(mu))) %*% exp(d %*% theta + x %*% beta)%*% d}) %>%
    reduce(`+`)/nrow(X)
}

psi_plpr_op <- function(theta, Y, D, X, mu, beta){
  y_list <- split(Y, seq(nrow(Y)))
  x_list <- split(X, seq(nrow(X)))
  d_list <- split(D, seq(nrow(D)))

  list(y_list, x_list, d_list) %>%
    pmap(function(y, x, d) {
      op <- t(d - x %*% t(mu))%*% (y - exp(d %*% theta + x %*% beta))
      op %*% t(op)}) %>%
    reduce(`+`)/nrow(Y)
}
