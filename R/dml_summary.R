#' @importFrom magrittr %>%
#' @importFrom tibble tibble as_tibble enframe
#' @importFrom broom glance tidy
#' @importFrom stats coef vcov fitted
#' @importFrom dplyr mutate
#'
#'
#'@method coef dml
#' @export
coef.dml <- function(x) x$coefficients

#'@method vcov dml
#' @export
vcov.dml <- function(x) x$vcov

#'@method fitted dml
#' @export
fitted.dml <- function(x) rep(0, x$nobs)

#'@method glance dml
#' @export
glance.dml <- function(x) tibble("n" = x$nobs,
                                 "r.squared" = NA_real_,
                                 "adj.r.squared" = NA_real_)

#'@method print dml
#' @export
print.dml <- function (x, digits = max(3L, getOption("digits") - 3L),
                       ...){
  if (length(coef(x))) {
    print.default(format(coef(x), digits = digits, big.mark = ","),
                  print.gap = 2L, quote = FALSE)
  }
  else cat("No coefficients\n")

  cat("\n")
  invisible(x)
}

#'@method summary dml
#' @export
summary.dml <- function(object){
  dml_sum <- list()
  dml_sum$call <- object$call
  dml_sum$coefficients <-
    cbind(Estimate = coef(object),
          `Std. Error` = sqrt(diag(vcov(object))))
  dml_sum$n <- object$nobs

  class(dml_sum) <- "summary.dml"
  dml_sum
}

#' @method tidy dml
#' @export
tidy.dml <- function(object, ...){
  tibble(term = names(object$coefficients),
         estimate = coef(object),
         std.error = sqrt(diag(vcov(object))))
}

#' @method print summary.dml
#' @export
print.summary.dml <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")

  cat("\nCoefficients:\n")
  if(length(x$coefficients) == 0)
    cat("(No coefficients)\n")
  else {
    printCoefmat(x$coefficients, digits = digits, big.mark = ",")
    cat("\nNumber of Observations:", format(x$n, big.mark = ","))
    cat("\n\n")
  }
}

