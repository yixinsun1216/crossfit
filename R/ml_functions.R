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
regression_forest2 <- function(f, d, num.trees = 1000, honesty = FALSE,
                               honesty.fraction = NULL, tune.parameters = TRUE, ...) {
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

  ff <- regression_forest(X, Y, num.trees = num.trees, honesty = honesty,
                          honesty.fraction = honesty.fraction,
                          tune.parameters = tune.parameters,...)

  ff[["formula"]] <- f
  class(ff) <- c("regression_forest", "grf")

  return(ff)
}

#' @export
#' @rdname predict_rlasso2
predict_rlasso2 <- function(rl, newdata = NULL) {
  f <- rl[["model"]]

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
    paste(y, "~ poly(", ., ", degree = ", deg, ")") %>%
    as.formula
  return(f_poly)
}

#' @export
run_ml <- function(f_gamma, fs_delta, folds, ml_type, dml_seed,
                   poly_degree, ...){
  # specify which ml method to use for training and predicting
  if(ml_type == "regression_forest"){
    train_ml <- "regression_forest2"
    predict_ml <- "predict_rf2"
  }

  if(ml_type == "lasso"){
    train_ml <- "cv.glmnet"
    predict_ml <- "predict"
  }

  if(ml_type == "rlasso"){
    train_ml <- "rlasso"
    predict_ml <- "predict"
  }

  if(ml_type %in% c("rlasso", "lasso")){
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

  set.seed(dml_seed)

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

  params <- get.params(...)
  print(paste("run_ml parameters", params))

  return(list(gamma = gamma, delta = delta))
}

