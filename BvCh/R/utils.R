## Internal helpers -------------------------------------------------------

#' Recycle arguments to a common length
#' @param ... Numeric or NULL arguments.
#' @return A list of arguments of equal (non-zero) length.
#' @keywords internal
#' @noRd
.recycle <- function(...) {
  args <- list(...)
  n <- 0L
  for (a in args) if (!is.null(a)) n <- max(n, length(a))
  lapply(args, function(a) {
    if (is.null(a)) return(NULL)
    if (length(a) == n) return(a)
    if (length(a) == 1L) return(rep(a, n))
    stop("All vector arguments must have length 1 or a common length.",
         call. = FALSE)
  })
}

#' Check that all parameters are finite and strictly positive
#' @keywords internal
#' @noRd
.chk_pos <- function(..., tol = .Machine$double.eps^0.5) {
  a <- c(...)
  if (!all(is.finite(a))) stop("Parameters must be finite.", call. = FALSE)
  if (any(a <= 0)) stop("Parameters must be strictly positive.", call. = FALSE)
  invisible(a)
}

#' Exact equality branch for two times, tolerant to floating point noise
#'
#' `z1 == z2` is the singular support of Marshall-Olkin type models. Data
#' read from files may differ in the last bits, so values within a relative
#' tolerance are treated as ties.
#' @keywords internal
#' @noRd
.is_tie <- function(z1, z2, tol = .Machine$double.eps^0.5) {
  abs(z1 - z2) <= tol * pmax(abs(z1), abs(z2))
}

#' Univariate Chen density (reliaR backend)
#'
#' Thin wrapper mapping the bivariate papers' scale parameter \eqn{\alpha}
#' onto the \code{lambda} argument of \code{\link[reliaR]{dchen}}.
#' @keywords internal
#' @noRd
.dChen <- function(x, beta, alpha) reliaR::dchen(x, beta, alpha)

#' Univariate Chen survival function (reliaR backend)
#' @keywords internal
#' @noRd
.sChen <- function(x, beta, alpha) reliaR::schen(x, beta, alpha)

#' Univariate Chen hazard function (reliaR backend)
#' @keywords internal
#' @noRd
.hChen <- function(x, beta, alpha) reliaR::hchen(x, beta, alpha)

#' @keywords internal
#' @noRd
.pChen <- function(q, beta, alpha) {
  ## reliaR::pchen rejects q <= 0, but F(0) = 0 for a Chen distribution, so
  ## clamp to the smallest positive double where q^beta underflows to 0 and
  ## the Chen CDF evaluates to exactly 0.
  reliaR::pchen(pmax(q, .Machine$double.xmin), beta, alpha)
}

#' Univariate Chen distribution function (reliaR backend)
#' @keywords internal
#' @noRd
.pChen_raw <- function(q, beta, alpha) reliaR::pchen(q, beta, alpha)

#' Univariate Chen quantile function (reliaR backend)
#' @keywords internal
#' @noRd
.qChen <- function(p, beta, alpha) reliaR::qchen(p, beta, alpha)

#' Univariate Chen random generation (reliaR backend)
#' @keywords internal
#' @noRd
.rChen <- function(n, beta, alpha) reliaR::qchen(stats::runif(n), beta, alpha)

#' log(1 + lambda * (expm1(y))) computed stably
#' @keywords internal
#' @noRd
.log1p_lam_expm1 <- function(y, lambda) {
  n <- max(length(y), length(lambda))
  y <- rep_len(y, n)
  lambda <- rep_len(lambda, n)
  v <- lambda * expm1(y)
  ## log1p(v) for v > 0 (the extended Chen support), v otherwise
  ifelse(v > 0, log1p(v), v)
}
