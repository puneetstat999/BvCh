#' Univariate Extended Chen distribution
#'
#' Density, distribution function, quantile function, random generation,
#' survival and hazard function for the univariate Extended Chen distribution
#' used as the latent building block of the Bivariate Extended Chen (BEC)
#' model.
#'
#' @details The Extended Chen survival function is
#' \deqn{S_{EC}(x;\alpha,\beta,\lambda) =
#'   \left[1+\lambda\left(\exp(x^{\beta})-1\right)\right]^{-\alpha},}
#' giving hazard
#' \deqn{h_{EC}(x;\alpha,\beta,\lambda) =
#'   \frac{\alpha\lambda\beta x^{\beta-1}\exp(x^{\beta})}
#'        {1+\lambda\left(\exp(x^{\beta})-1\right)}}
#' and density \eqn{f_{EC}(x) = S_{EC}(x)\,h_{EC}(x)}. The quantile function is
#' available in closed form,
#' \deqn{Q(p) = \left[\log\left\{1+\frac{p^{-1/\alpha}-1}{\lambda}\right\}\right]^{1/\beta}.}
#'
#' Unlike the univariate Chen distribution, which \pkg{reliaR} supplies, no
#' Extended Chen implementation exists there, so it is provided here.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param alpha,lambda,beta positive parameters.
#' @param log,log.p logical; if \code{TRUE} probabilities are returned as
#'   \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default) \code{P[X <= x]},
#'   otherwise \code{P[X > x]}.
#'
#' @return \code{dech} gives the density, \code{pech} the distribution
#'   function, \code{qech} the quantile function, \code{rech} random
#'   deviates, \code{sech} the survival function and \code{hech} the hazard
#'   function.
#'
#' @examples
#' sech(0.5, alpha = 1.5, beta = 1.5, lambda = 2)
#' hech(0.5, alpha = 1.5, beta = 1.5, lambda = 2)
#' qech(c(0.25, 0.5, 0.75), alpha = 1.5, beta = 1.5, lambda = 2)
#'
#' @name extendedchen
NULL

#' @rdname extendedchen
#' @export
sech <- function(x, alpha, beta, lambda) {
  .chk_pos(alpha, beta, lambda)
  pmin((1 + lambda * expm1(pmax(x, 0)^beta))^(-alpha), 1)
}

#' @rdname extendedchen
#' @export
hech <- function(x, alpha, beta, lambda) {
  .chk_pos(alpha, beta, lambda)
  x <- pmax(x, 0)
  y <- x^beta
  (alpha * lambda * beta * x^(beta - 1) * exp(y)) / (1 + lambda * expm1(y))
}

#' @rdname extendedchen
#' @export
dech <- function(x, alpha, beta, lambda, log = FALSE) {
  f <- sech(x, alpha, beta, lambda) * hech(x, alpha, beta, lambda)
  f[!is.finite(f)] <- 0
  if (log) log(f) else f
}

#' @rdname extendedchen
#' @export
pech <- function(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
  s <- sech(q, alpha, beta, lambda)
  p <- if (lower.tail) 1 - s else s
  if (log.p) log(p) else p
}

#' @rdname extendedchen
#' @export
qech <- function(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
  .chk_pos(alpha, beta, lambda)
  if (log.p) p <- exp(p)
  if (any(p < 0 | p > 1)) stop("Probabilities must lie in [0, 1].", call. = FALSE)
  ## s is the survival probability we need to invert
  s <- if (lower.tail) 1 - p else p
  out <- rep(0, length(s))
  ok <- s > 0
  out[ok] <- log1p((s[ok]^(-1 / alpha) - 1) / lambda)^(1 / beta)
  out
}

#' @rdname extendedchen
#' @export
rech <- function(n, alpha, beta, lambda) {
  .chk_pos(alpha, beta, lambda)
  qech(stats::runif(n), alpha, beta, lambda)
}
