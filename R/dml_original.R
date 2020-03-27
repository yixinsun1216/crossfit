dml_step_original <- function(f, d, model, dml_seed = NULL,
                              ml, poly_degree, ...){
   # ensure fm has no intercept in the parametric part
  f <-
    Formula(f) %>%
    update(. ~ 0 + . | 0 + .)

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
  folds$train <- map(folds$train, data.frame)

  # 1. Calculate theta and beta using training data
  # 2. Use training data and estimates from 1 to find mu
  # 3. Optimize over psi using hold out data
  f_ml <- paste(c(names(tx), names(td)), collapse = " + ") %>%
    paste0(names(ty), " ~ ", .) %>%
    as.formula

  # Step 1:
  ml_coefs <- estimate_ml(f_ml, folds, ml, poly_degree, ...)

  # Step 2: Use estimated coefficients to find
  # mu = E[X'exp(D*theta + X*beta)X]^{-1}E[X'exp(D*theta - X*beta)D]
  mu_hat <-
    map2(ml_coefs, folds$train,
          function(x, y) estimate_mu(x, y, names(td), names(tx)))

  # Step 3:

}


# ============================================================================
# Use Original DML approach (not the concentrating out approach)
# ============================================================================
estimate_ml <- function(f_ml, folds, ml_type, poly_degree, ...){
  if(ml_type == "regression_forest"){
    train_ml <- "regression_forest2"
  }

  if(ml_type == "lasso"){
    train_ml <- "cv.glmnet"
  }

  if(ml_type == "rlasso"){
    train_ml <- "rlasso2"
  }

  if(ml_type %in% c("rlasso", "lasso") & as.numeric(poly_degree) > 1){
    # create new formulas that holds all the permutation of polynomials
    # basis functions for lasso
    f_ml <- poly_formula(f_ml, poly_degree)
  }

  # train and estimate coefficients beta_hat and theta_hat for the formula
  # Y = X*Beta_hat + D*theta_hat
  map(folds$train, function(y) {
      ml_coef <- coef(get(train_ml)(f_ml, y, ...))
      data.frame(coef = pluck(dimnames(ml_coef), 1), value = matrix(ml_coef)) %>%
        mutate(coef = as.character(coef)) %>%
        mutate(coef = if_else(str_detect(coef, regex("Intercept", ignore.case = TRUE)), "Intercept", coef))
  })
}

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

  return(solve(j_bb) %*% j_tb)
}

j_calc <- function(X, D, theta, beta, int, type){
  V <- D
  if(type != "tb") V <- X
  t(as.matrix(X*exp(sum(D[1,]*theta) + sum(X*beta) + int))) %*% as.matrix(V)
}
