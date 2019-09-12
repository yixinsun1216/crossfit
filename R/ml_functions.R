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
#' @importFrom Formula Formula
#'
#'

default_forest <- list(num.trees = 1000,
                       honesty = TRUE,
                       honesty.fraction = NULL)

#' @export
regression_forest2 <- function(f, d, num.trees = 1000, honesty = TRUE,
                               honesty.fraction = NULL) {
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

  ff <- regression_forest(X, Y, default_forest)

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

cv.glmnet2 <- function(f, d, ...){

  lf <- cv.glmnet(f, data = d)

}
