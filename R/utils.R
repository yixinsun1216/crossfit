#' Auxilary functions to main dml process
#'
#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble enframe
#' @importFrom purrr pluck map2
#' @importFrom stats update formula model.matrix model.frame predict
#' @importFrom grf regression_forest
#' @importFrom Formula Formula


# allows formulas to be passsed to grf::regression_forest instead of matrices
regression_forest2 <- function(f, d, args) {
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

  rf_args <- append(list(X = X, Y = Y, num.trees = 1000), args)
  rf_args <- rf_args[!duplicated(names(rf_args), fromLast = TRUE)]
  ff <- do.call(regression_forest, rf_args)

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

square <- function(x) 0.5 * t(x) %*% x


get_rhs_cols <- function(f, d, part = 1) {
  options(na.action='na.pass')

  f %>%
    formula(lhs = 0, rhs = part) %>%
    model.matrix(d) %>%
    as_tibble() %>%
    rename_all(~ str_replace_all(., regex("[(, =)]"), "_"))
}

get_lhs_col <- function(f, d) {
  f %>%
    formula(lhs = 1, rhs = 0) %>%
    model.frame(d, na.action = NULL) %>%
    as_tibble %>%
    rename_all(~ str_replace_all(., regex("[(, =)]"), "_"))
}


