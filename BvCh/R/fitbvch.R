#' Maximum likelihood estimation for the BvCh model
#'
#' Fits the Marshall-Olkin type Bivariate Chen model by maximum likelihood.
#'
#' @details Observations with \code{z1 == z2} fall on the singular component
#'   of the model. Whether they contribute the singular density
#'   \eqn{f_0} or are discarded is controlled by \code{singular}; the counts
#'   are always reported so the choice stays visible. Standard errors come
#'   from the numerical Hessian of the log-likelihood and are on the
#'   original (untransformed) scale.
#'
#' @param x1,x2 numeric vectors of observations, or a two column matrix.
#' @param start optional numeric vector of starting values on the natural
#'   scale, in the order \code{alpha1, alpha2, alpha3, beta}.
#' @param singular logical; whether tied observations contribute the
#'   singular density.
#' @param method optimisation method passed to \code{\link[stats]{optim}}.
#' @param control optional control list for \code{\link[stats]{optim}}.
#' @param se logical; whether standard errors should be computed.
#'
#' @return An object of class \code{fitbvch}: a list with components
#'   \code{coefficients}, \code{vcov}, \code{se}, \code{loglik},
#'   \code{aic}, \code{bic}, \code{nobs}, \code{counts}, \code{converged},
#'   \code{message} and \code{call}.
#'
#' @seealso \code{\link{fitbec}}, \code{\link{fitifgmchen}},
#'   \code{\link{logLik_bvch}}
#'
#' @examples
#' set.seed(1)
#' d <- rbvch(400, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#' fit <- fitbvch(d)
#' fit
#' coef(fit)
#'
#' @export
fitbvch <- function(x1, x2 = NULL, start = NULL, singular = TRUE,
                    method = "BFGS", control = list(), se = TRUE) {
  cl <- match.call()
  d <- .as_pair(x1, x2, "fitbvch")
  z1 <- d[, 1]; z2 <- d[, 2]
  st <- if (is.null(start)) .default_start("bvch") else log(start)
  if (length(st) != 4L) stop("'start' must have length 4.", call. = FALSE)
  tie <- .is_tie(z1, z2)
  use <- if (singular) rep(TRUE, length(z1)) else !tie
  if (!any(use)) stop("No usable observations.", call. = FALSE)
  f <- function(p) .nll(p, z1[use], z2[use], "bvch", singular)
  op <- stats::optim(st, f, method = method, control = control)
  est <- exp(op$par)
  names(est) <- c("alpha1", "alpha2", "alpha3", "beta")
  ll <- -op$value
  ab <- .aic_bic(ll, sum(use), 4L)
  V <- if (se && op$convergence == 0)
    .hessian_vcov(op$par, f, TRUE, 4L) else matrix(NA_real_, 4, 4)
  structure(list(
    coefficients = est, vcov = V, se = sqrt(pmax(diag(V), 0)),
    loglik = ll, aic = ab$aic, bic = ab$bic, negloglik = ab$negloglik,
    nobs = sum(use),
    counts = c(.mo_counts(z1, z2, tie), used = sum(use)),
    converged = op$convergence == 0,
    message = op$message, model = "bvch", call = cl),
    class = "fitbvch")
}
