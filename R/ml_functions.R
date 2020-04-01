#' Machine Learning Functions
#'
#' \code{regression_forest2} estimates a regression forest with a formula
#'   instead of manually passing matrices.
#' \code{predict_rf2} predicts based on the forest ouputted from
#'   \code{regression_forest2}.

#' @param f an object of class formula representing the model to be fitted.
#' @param d a dataframe containing the variables in f.
#' @param ... additional arguments taht can be passed to
#'   \code{\link[grf:regression_forest]{regression_forest}},such as num.trees,
#'   honesty,honesty.fraction,seed.
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble enframe
#' @importFrom purrr pluck map2
#' @importFrom stats update formula model.matrix model.frame predict
#' @importFrom grf regression_forest
#' @importFrom glmnetUtils cv.glmnet
#' @importFrom hdm rlasso
#' @importFrom Formula Formula
#'
#'

#' @export
# note on passing in arguments https://stackoverflow.com/questions/23922130/default-argument-in-r-function-formal-argument-matched-by-multiple-actual-argum?rq=1
regression_forest2 <- function(f, d, num.trees = 1000, ...) {
  f <- Formula(f)

  Y <-
    f %>%
    formula(lhs = 1, rhs = 0) %>%
    model.frame(d) %>%
    as.matrix

  X <-
    formula(f, rhs = 1, lhs = 0) %>%
    update(~ 0 + .) %>%
    model.matrix(d)

  ff <- regression_forest(X, Y, num.trees = num.trees,...)

  ff[["formula"]] <- f
  class(ff) <- c("regression_forest", "grf")

  return(ff)
}

#' @export
# wrapper to return the original formula as a part of the rlasso object
rlasso2 <- function(f, d, ...){
  rlasso_model <- rlasso(f, d, ...)
  rlasso_model[["formula"]] <- f
 return(rlasso_model)
}

#' @export
predict_rlasso2 <- function(rl, newdata = NULL) {
  f_rl <- rl[["formula"]]

  if(!is.null(newdata)) {
    X <-
      formula(f_rl, rhs = 1, lhs = 0) %>%
      update(~ 0 + .) %>%
      model.matrix(newdata)
    return(predict(rl, X))

  } else {
    return(predict(rl))
  }
}

predict_rf2 <- function(forest, newdata = NULL) {
  f <- forest[["formula"]]

  if(!is.null(newdata)) {
    X <-
      formula(f, rhs = 1, lhs = 0) %>%
      update(~ 0 + .) %>%
      model.matrix(newdata)
    return(pluck(predict(forest, X), "predictions"))
  } else {
    return(pluck(predict(forest), "predictions"))
  }
}

# expand a formula to all polynomial terms
poly_formula <- function(f, deg){
  y <- all.vars(f)[1]
  f_poly <-
    labels(terms(f)) %>%
    paste(collapse = " , ") %>%
    paste(y, "~ poly(", ., ", raw = FALSE, degree = ", deg, ")") %>%
    as.formula
  return(f_poly)
}

#' @export
run_ml <- function(f_gamma, fs_delta, folds, ml, poly_degree, ...){
  # specify which ml method to use for training and predicting
  if(ml == "regression_forest"){
    train_ml <- "regression_forest2"
    predict_ml <- "predict_rf2"
  }

  if(ml == "lasso"){
    train_ml <- "cv.glmnet"
    predict_ml <- "predict"
  }

  if(ml == "rlasso"){
    train_ml <- "rlasso2"
    predict_ml <- "predict_rlasso2"
  }

  if(ml %in% c("rlasso", "lasso")){
    # create new formulas that holds all the permutation of polynomials
    # basis functions for lasso
    if(as.numeric(poly_degree) > 0){
      f_gamma <- poly_formula(f_gamma, poly_degree)
      fs_delta <- map(fs_delta, poly_formula, poly_degree)
    }

    # shape dataframe so it is usable with glmnet
    folds$test <- map(folds$test, data.frame)
    folds$train <- map(folds$train, data.frame)
  }

  # train and estimate values of delta in hold out sample
  gamma_models <- map(folds$train, function(y) {
    get(train_ml)(f_gamma, y, ...)
  })

  gamma <-
    map2(gamma_models, folds$test, ~ get(predict_ml)(.x, .y)) %>%
    map(as.matrix)

  # train and estimate values of delta in the hold out sample
  delta_models <-
    folds$train %>%
    map(function(x) map(fs_delta, function(y) {
      get(train_ml)(y, x, ...)
    }))

  delta <-
    map2(delta_models, folds$test,
         function(x, y)
           map2(x, list(y), ~ get(predict_ml)(.x, .y)) %>%
           unlist %>%
           matrix(ncol = length(fs_delta)))

  return(list(gamma = gamma, delta = delta))
}

