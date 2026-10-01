#' Maximum likelihood estimation for the BEC model
#'
#' Fits the Bivariate Extended Chen model by maximum likelihood.
#'
#' @inheritParams fitbvch
#' @param start optional starting values on the natural scale in the order
#'   \code{alpha1, alpha2, alpha3, beta, lambda}.
#'
#' @return An object of class \code{fitbec} with the same components as
#'   \code{\link{fitbvch}}.
#'
#' @seealso \code{\link{fitbvch}}, \code{\link{logLik_bec}}
#'
#' @examples
#' set.seed(1)
#' d <- rbec(400, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8,
#'           lambda = 2)
#' fit <- fitbec(d)
#' fit
#'
#' @export
fitbec <- function(x1, x2 = NULL, start = NULL, singular = TRUE,
                   method = "BFGS", control = list(), se = TRUE) {
  cl <- match.call()
  d <- .as_pair(x1, x2, "fitbec")
  z1 <- d[, 1]; z2 <- d[, 2]
  st <- if (is.null(start)) .default_start("bec") else log(start)
  if (length(st) != 5L) stop("'start' must have length 5.", call. = FALSE)
  tie <- .is_tie(z1, z2)
  use <- if (singular) rep(TRUE, length(z1)) else !tie
  if (!any(use)) stop("No usable observations.", call. = FALSE)
  f <- function(p) .nll(p, z1[use], z2[use], "bec", singular)
  op <- stats::optim(st, f, method = method, control = control)
  est <- exp(op$par)
  names(est) <- c("alpha1", "alpha2", "alpha3", "beta", "lambda")
  ll <- -op$value
  ab <- .aic_bic(ll, sum(use), 5L)
  V <- if (se && op$convergence == 0)
    .hessian_vcov(op$par, f, TRUE, 5L) else matrix(NA_real_, 5, 5)
  structure(list(
    coefficients = est, vcov = V, se = sqrt(pmax(diag(V), 0)),
    loglik = ll, aic = ab$aic, bic = ab$bic, negloglik = ab$negloglik,
    nobs = sum(use),
    counts = c(.mo_counts(z1, z2, tie), used = sum(use)),
    converged = op$convergence == 0, message = op$message,
    model = "bec", call = cl),
    class = "fitbec")
}
