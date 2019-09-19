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
#' @importFrom purrr pluck
#' @importFrom stats update formula model.matrix model.frame predict
#' @importFrom grf regression_forest
#' @importFrom glmnetUtils cv.glmnet
#' @importFrom Formula Formula
#'
#'

#' @export
# note on passing in arguments https://stackoverflow.com/questions/23922130/default-argument-in-r-function-formal-argument-matched-by-multiple-actual-argum?rq=1
regression_forest2 <- function(f, d, num.trees = 1000, honesty = TRUE,
                               honesty.fraction = NULL, ...) {
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
                          honesty.fraction = honesty.fraction, ...)

  ff[["formula"]] <- f
  class(ff) <- c("regression_forest", "grf")

  return(ff)
}

#' @export
#' @rdname regression_forest2
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



#' @export
run_ml <- function(f, folds_test, ml_type, value_type, dml_seed, ...){
  # specify which ml method to use for training and predicting
  if(ml_type == "regression_forest"){
    train_ml <- "regression_forest2"
    predict_ml <- "predict_rf2"
  }

  if(ml_type == "lasso"){
    train_ml <- "cv.glmnet"
    predict_ml <- "predict"
    folds_test <- map(folds_test, data.frame)
  }

  set.seed(dml_seed)

  # train and estimate values of delta in hold out sample
  if(value_type == "gamma"){
    gamma_models <- map(folds_test, function(y) {
      get(train_ml)(f, d = y, ...)
    })

    output <-
      map2(gamma_models, folds_test, ~ get(predict_ml)(.x, .y)) %>%
      map(as.matrix)

  }

  # train and estimate values of gamma in the hold out sample
  if(value_type == "delta"){
    delta_models <-
      folds_test %>%
      map(function(x) map(f, function(y) {
        get(train_ml)(y, d = x, ...)
      }))

    output <-
      map2(delta_models, folds_test,
           function(x, y)
             map2(x, list(y), ~ get(predict_ml)(.x, .y)) %>%
             unlist %>%
             matrix(ncol = length(f)))
  }

  return(output)
}
