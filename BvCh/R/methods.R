## S3 methods for the fitted model objects.  Each method block carries
## `@rdname print.fitbvch` and `@export`, so roxygen2 registers the
## corresponding S3method() directives; the help topic itself is the
## documentation block at the end of this file.

#' @rdname print.fitbvch
#' @export
print.fitbvch <- function(x, digits = 4, ...) .print_fit(x, digits)

#' @rdname print.fitbvch
#' @export
print.fitbec <- function(x, digits = 4, ...) .print_fit(x, digits)

#' @rdname print.fitbvch
#' @export
print.fitifgmchen <- function(x, digits = 4, ...) .print_fit(x, digits)

#' @rdname print.fitbvch
#' @export
summary.fitbvch <- function(object, digits = 4, signif.stars = FALSE, ...) {
  .print_fit(object, digits)
  invisible(object)
}

#' @rdname print.fitbvch
#' @export
summary.fitbec <- summary.fitbvch

#' @rdname print.fitbvch
#' @export
summary.fitifgmchen <- summary.fitbvch

#' @rdname print.fitbvch
#' @export
coef.fitbvch <- function(object, ...) object$coefficients

#' @rdname print.fitbvch
#' @export
coef.fitbec <- coef.fitbvch

#' @rdname print.fitbvch
#' @export
coef.fitifgmchen <- coef.fitbvch

#' @rdname print.fitbvch
#' @export
vcov.fitbvch <- function(object, ...) object$vcov

#' @rdname print.fitbvch
#' @export
vcov.fitbec <- vcov.fitbvch

#' @rdname print.fitbvch
#' @export
vcov.fitifgmchen <- vcov.fitbvch

#' @rdname print.fitbvch
#' @export
nobs.fitbvch <- function(object, ...) object$nobs

#' @rdname print.fitbvch
#' @export
nobs.fitbec <- nobs.fitbvch

#' @rdname print.fitbvch
#' @export
nobs.fitifgmchen <- nobs.fitbvch

#' @rdname print.fitbvch
#' @export
logLik.fitbvch <- function(object, ...) {
  structure(object$loglik, df = 4L, nobs = object$nobs)
}

#' @rdname print.fitbvch
#' @export
logLik.fitbec <- function(object, ...) {
  structure(object$loglik, df = 5L, nobs = object$nobs)
}

#' @rdname print.fitbvch
#' @export
logLik.fitifgmchen <- function(object, ...) {
  structure(object$loglik, df = 6L, nobs = object$nobs)
}

.print_fit <- function(x, digits = 4) {
  lab <- switch(x$model,
    bvch = "BvCh", bec = "BEC", ifgmchen = "IFGM-Chen", x$model)
  cat("\n", lab, " model fit by maximum likelihood\n\n", sep = "")
  cat("Call: ", paste(deparse(x$call), collapse = " "), "\n", sep = "")
  cat("\nCoefficients:\n")
  out <- data.frame(Estimate = x$coefficients, `Std. Error` = x$se,
                    check.names = FALSE)
  print(round(out, digits))
  cat("\nn = ", x$nobs,
      "   (z1 < z2: ", x$counts[["below"]],
      ", z1 > z2: ", x$counts[["above"]],
      ", z1 = z2: ", x$counts[["tie"]], ")\n", sep = "")
  cat("-log L = ", format(x$negloglik, digits = digits),
      "   AIC = ", format(x$aic, digits = digits),
      "   BIC = ", format(x$bic, digits = digits), "\n", sep = "")
  if (!x$converged)
    cat("Warning: optimiser did not converge:", x$message, "\n")
  invisible(x)
}

#' Print, summary and extractor methods for fitted bivariate Chen models
#'
#' Methods for the objects returned by \code{\link{fitbvch}},
#' \code{\link{fitbec}} and \code{\link{fitifgmchen}}.
#'
#' @details The fitted objects also carry \code{aic} and \code{bic}
#'   components, and \code{\link{compare_models}} tabulates them across
#'   models. Because \code{aic} and \code{bic} are not S3 generics that a
#'   package can register methods for, no \code{aic} or \code{bic} methods
#'   are provided.
#'
#' @param x,object a fitted model object.
#' @param digits number of significant digits to print.
#' @param signif.stars logical; whether to print significance stars.
#' @param ... further arguments, currently unused.
#'
#' @return The printing methods return their argument invisibly;
#'   \code{coef}, \code{vcov}, \code{nobs} and \code{logLik} return the
#'   corresponding component of the fit.
#'
#' @usage
#' \method{print}{fitbvch}(x, digits = 4, ...)
#' \method{print}{fitbec}(x, digits = 4, ...)
#' \method{print}{fitifgmchen}(x, digits = 4, ...)
#' \method{summary}{fitbvch}(object, digits = 4, signif.stars = FALSE, ...)
#' \method{coef}{fitbvch}(object, ...)
#' \method{coef}{fitbec}(object, ...)
#' \method{coef}{fitifgmchen}(object, ...)
#' \method{vcov}{fitbvch}(object, ...)
#' \method{vcov}{fitbec}(object, ...)
#' \method{vcov}{fitifgmchen}(object, ...)
#' \method{nobs}{fitbvch}(object, ...)
#' \method{nobs}{fitbec}(object, ...)
#' \method{nobs}{fitifgmchen}(object, ...)
#' \method{logLik}{fitbvch}(object, ...)
#' \method{logLik}{fitbec}(object, ...)
#' \method{logLik}{fitifgmchen}(object, ...)
#'
#' @examples
#' set.seed(1)
#' d <- rbvch(300, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#' fit <- fitbvch(d)
#' print(fit)
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#' nobs(fit)
#'
#' @name print.fitbvch
NULL
