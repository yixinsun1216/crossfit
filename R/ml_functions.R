
# estimate a regression forest with a formula instead of manually passing
# matrices.  the formula y ~ x1 + x2 + x3 means estimate a regression forest
# for the outcome y, where the splitting variables are x1, x2, and x3
# to get predicted values, pass the results of this to predict()
regression_forest2 <- function(f, d, ...) {
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

  ff <- regression_forest(X, Y, ...)

  ff[["formula"]] <- f
  class(ff) <- c("regression_forest", "grf")

  return(ff)
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
