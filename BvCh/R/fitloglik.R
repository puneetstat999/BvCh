#' Log-likelihood functions of the bivariate Chen models
#'
#' Log-likelihoods used by \code{\link{fitbvch}}, \code{\link{fitbec}} and
#' \code{\link{fitifgmchen}}.
#'
#' @details For the Marshall-Olkin type models the log-likelihood is the
#'   exact observed-data log-likelihood: the contribution of an observation
#'   is entirely determined by whether \eqn{z_1<z_2}, \eqn{z_1>z_2} or
#'   \eqn{z_1=z_2}, so the complete-data log-likelihood coincides with the
#'   observed-data one. Direct numerical maximisation is therefore used;
#'   it reaches the same fixed point as the EM algorithm.
#'
#' @param z1,z2 numeric vectors of observations.
#' @param par numeric vector of parameters on the log scale.
#' @param singular logical; if \code{TRUE} ties (\code{z1 == z2}) contribute
#'   the singular density. Use \code{FALSE} to drop tied observations.
#'
#' @return The sum of the log densities, or \code{-Inf} when any contribution
#'   is not finite.
#'
#' @examples
#' set.seed(1)
#' d <- rbvch(50, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#' logLik_bvch(log(c(1.5, 2.5, 0.5, 1.8)), d[, 1], d[, 2])
#'
#' @name loglik
NULL

#' @keywords internal
#' @noRd
.loglik_impl <- function(par, z1, z2, model, singular = TRUE) {
  ## model is one of "bvch", "bec" or "ifgmchen"
  p <- exp(par)
  ## the optimiser may propose extreme values; reject them with a bad
  ## objective rather than an error so the search can continue
  if (any(!is.finite(p)) || any(p <= 0)) return(-Inf)
  if (model == "bvch") {
    ll <- dbvch(z1, z2, p[1], p[2], p[3], p[4], log = TRUE)
  } else if (model == "bec") {
    ll <- dbec(z1, z2, p[1], p[2], p[3], p[4], p[5], log = TRUE)
  } else {
    ll <- difgmchen(z1, z2, p[1], p[2], p[3], p[4], p[5], p[6], log = TRUE)
  }
  if (!singular) ll <- ll[!.is_tie(z1, z2)]
  if (any(!is.finite(ll))) return(-Inf)
  sum(ll)
}

#' @rdname loglik
#' @export
logLik_bvch <- function(par, z1, z2, singular = TRUE) {
  .loglik_impl(par, z1, z2, "bvch", singular)
}

#' @rdname loglik
#' @export
logLik_bec <- function(par, z1, z2, singular = TRUE) {
  .loglik_impl(par, z1, z2, "bec", singular)
}

#' @rdname loglik
#' @export
logLik_ifgmchen <- function(par, z1, z2) {
  .loglik_impl(par, z1, z2, "ifgmchen", TRUE)
}

#' @keywords internal
#' @noRd
.aic_bic <- function(loglik, n, k) {
  list(aic = 2 * k - 2 * loglik, bic = k * log(n) - 2 * loglik,
       negloglik = -loglik)
}

#' @keywords internal
#' @noRd
.nll <- function(par, z1, z2, model, singular) {
  v <- .loglik_impl(par, z1, z2, model, singular)
  if (!is.finite(v)) 1e300 else -v
}

#' @keywords internal
#' @noRd
.default_start <- function(model) {
  switch(model,
    bvch = log(c(1.0, 1.0, 1.0, 1.5)),
    bec = log(c(1.0, 1.0, 1.0, 1.5, 1.0)),
    ifgmchen = log(c(1.0, 1.0, 1.0, 1.0, 0.5, 0.5)),
    log(c(1, 1, 1, 1)))
}

#' @keywords internal
#' @noRd
.as_pair <- function(x1, x2 = NULL, fn = "fit") {
  if (is.null(x2) && is.matrix(x1) && ncol(x1) == 2L) {
    m <- x1
  } else if (is.null(x2)) {
    stop(paste0(fn, ": supply either a two column matrix or two vectors."),
         call. = FALSE)
  } else {
    if (length(x1) != length(x2))
      stop("'x1' and 'x2' must have the same length.", call. = FALSE)
    m <- cbind(x1, x2)
  }
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  if (any(!is.finite(m)))
    stop("Observations must be finite.", call. = FALSE)
  if (any(m < 0))
    stop("Observations must be non-negative.", call. = FALSE)
  if (nrow(m) < 2L)
    stop("At least two observations are required.", call. = FALSE)
  m
}

#' @keywords internal
#' @noRd
.hessian_vcov <- function(par, f, natural = TRUE, k = length(par)) {
  ## f is the NEGATIVE log-likelihood, so its Hessian is positive definite
  ## at the optimum and solving it gives the covariance matrix directly.
  H <- try(stats::optimHess(par, f), silent = TRUE)
  if (inherits(H, "try-error")) return(matrix(NA_real_, k, k))
  Vi <- try(solve(H), silent = TRUE)
  if (inherits(Vi, "try-error")) return(matrix(NA_real_, k, k))
  Vi <- as.matrix(Vi)
  if (natural) {
    ## delta method back to the natural scale
    diag(Vi) <- diag(Vi) * exp(par)^2
  }
  Vi
}

#' @keywords internal
#' @noRd
.mo_counts <- function(z1, z2, tie) {
  c(below = sum(z1 < z2), above = sum(z1 > z2), tie = sum(tie))
}
