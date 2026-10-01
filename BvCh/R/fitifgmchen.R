#' Maximum likelihood estimation for the IFGM-Chen model
#'
#' Fits the Bivariate Inverse FGM Chen model by maximum likelihood.
#'
#' @details The copula parameters \code{gamma} and \code{omega} are constrained
#'   to \eqn{[0,1]} and are estimated on the logit scale; the remaining
#'   parameters are estimated on the log scale. All reported values are on
#'   the natural scale.
#'
#' @section Identifiability:
#' The four marginal parameters are estimated very precisely, but
#'   \code{gamma} and \code{omega} are only weakly identified: with
#'   \eqn{n=5000} simulated from \eqn{(\gamma,\omega)=(0.6,0.4)} the marginal
#'   parameters have standard errors around \eqn{0.014} while \code{gamma}
#'   and \code{omega} have standard errors around \eqn{0.54} and
#'   \eqn{1.29}. The log-likelihood surface is nearly flat in that
#'   direction, so individual values of \code{gamma} and \code{omega} should
#'   not be over-interpreted; \code{logLik_ifgmchen} is exported so the
#'   profile likelihood can be inspected and \code{start} used to explore
#'   it.
#'
#' @inheritParams fitbvch
#' @param start optional starting values on the natural scale in the order
#'   \code{theta1, beta1, theta2, beta2, gamma, omega}.
#'
#' @return An object of class \code{fitifgmchen} with the same components as
#'   \code{\link{fitbvch}}.
#'
#' @seealso \code{\link{fitbvch}}, \code{\link{logLik_ifgmchen}}
#'
#' @examples
#' set.seed(1)
#' d <- rifgmchen(400, theta1 = 1.2, beta1 = 1.3, theta2 = 1.5,
#'                beta2 = 1.1, gamma = 0.6, omega = 0.4)
#' fit <- fitifgmchen(d)
#' fit
#'
#' @export
fitifgmchen <- function(x1, x2 = NULL, start = NULL, method = "BFGS",
                        control = list(), se = TRUE) {
  cl <- match.call()
  d <- .as_pair(x1, x2, "fitifgmchen")
  z1 <- d[, 1]; z2 <- d[, 2]
  st <- if (is.null(start)) .default_start("ifgmchen") else log(start)
  if (length(st) != 6L) stop("'start' must have length 6.", call. = FALSE)
  nll <- function(p) {
    q <- c(exp(p[1:4]), stats::plogis(p[5]), stats::plogis(p[6]))
    if (any(!is.finite(q)) || any(q <= 0)) return(1e300)
    v <- try(.loglik_impl(log(q), z1, z2, "ifgmchen", TRUE), silent = TRUE)
    if (inherits(v, "try-error") || !is.finite(v)) 1e300 else -v
  }
  q0 <- exp(st)
  st2 <- c(log(q0[1:4]), qlogis(q0[5]), qlogis(q0[6]))
  op <- stats::optim(st2, nll, method = method, control = control)
  est <- c(exp(op$par[1:4]), stats::plogis(op$par[5:6]))
  names(est) <- c("theta1", "beta1", "theta2", "beta2", "gamma", "omega")
  ll <- -op$value
  ab <- .aic_bic(ll, length(z1), 6L)
  V <- if (se && op$convergence == 0)
    .hessian_vcov(op$par, nll, FALSE, 6L) else matrix(NA_real_, 6, 6)
  structure(list(
    coefficients = est, vcov = V, se = sqrt(pmax(diag(V), 0)),
    loglik = ll, aic = ab$aic, bic = ab$bic, negloglik = ab$negloglik,
    nobs = length(z1),
    counts = c(.mo_counts(z1, z2, .is_tie(z1, z2)), used = length(z1)),
    converged = op$convergence == 0, message = op$message,
    model = "ifgmchen", call = cl),
    class = "fitifgmchen")
}
