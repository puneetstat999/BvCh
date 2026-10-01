#' Univariate Chen distribution
#'
#' Density, distribution function, quantile function, random generation,
#' survival and hazard function for the univariate Chen distribution with
#' shape parameter \code{beta} and scale parameter \code{lambda}.
#'
#' @details These are the \pkg{reliaR} implementations, re-exported by
#'   \pkg{BvCh} so that the Chen building blocks used in the bivariate
#'   formulae are available without attaching a second package. The density
#'   is
#'   \deqn{f(x;\lambda,\beta) = \lambda\beta x^{\beta-1}
#'     \exp\left(x^{\beta}\right)
#'     \exp\left[\lambda\left\{1-\exp\left(x^{\beta}\right)\right\}\right],}
#'   and the survival function is
#'   \deqn{S(x;\lambda,\beta) = \exp\left[\lambda\left\{1-\exp(x^{\beta})\right\}\right].}
#'   Taking \eqn{\lambda = \alpha} gives the \eqn{S_{Ch}(z;\alpha,\beta)} that
#'   appears in the bivariate Chen formulae.
#'
#' @section Random generation:
#' \code{rchen} is \emph{not} a direct call to \code{\link[reliaR]{rchen}}.
#' In \pkg{reliaR} 0.2 that function does not agree with the package's own
#' distribution, quantile and survival functions: for \eqn{\beta=1.5},
#' \eqn{\lambda=2} its output has Kolmogorov-Smirnov distance about 0.38 from
#' \code{pchen} (\eqn{p<10^{-16}}), and the discrepancy grows with
#' \eqn{\beta} (distance about 0.85 at \eqn{\beta=2.5}); it happens to be
#' correct only when \eqn{\beta=1}. Since \code{qchen} is the verified
#' analytic inverse of \code{pchen} (agreement to machine precision),
#' \code{rchen} here draws \code{qchen(runif(n), beta, lambda)}, which keeps
#' all six univariate Chen routines mutually consistent. Everything is still
#' computed by \pkg{reliaR}.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param beta shape parameter.
#' @param lambda scale parameter.
#' @param log,log.p logical; if \code{TRUE} probabilities are returned as
#'   \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default) \code{P[X <= x]},
#'   otherwise \code{P[X > x]}.
#' @param t age component for the conditional reliability function.
#'
#' @return \code{dchen} gives the density, \code{pchen} the distribution
#'   function, \code{qchen} the quantile function, \code{rchen} random
#'   deviates, \code{schen} the survival function and \code{hchen} the hazard
#'   function.
#'
#' @examples
#' dchen(0.5, beta = 1.5, lambda = 2)
#' pchen(0.5, beta = 1.5, lambda = 2)
#' qchen(0.5, beta = 1.5, lambda = 2)
#' schen(0.5, beta = 1.5, lambda = 2)
#' hchen(0.5, beta = 1.5, lambda = 2)
#'
#' @seealso \code{\link[reliaR]{Chen}}, \code{\link[reliaR]{Chensurvival}}
#' @name chen
NULL

#' @rdname chen
#' @export
dchen <- function(x, beta, lambda, log = FALSE) {
  reliaR::dchen(x, beta = beta, lambda = lambda, log = log)
}

#' @rdname chen
#' @export
pchen <- function(q, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
  reliaR::pchen(q, beta = beta, lambda = lambda,
                lower.tail = lower.tail, log.p = log.p)
}

#' @rdname chen
#' @export
qchen <- function(p, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
  reliaR::qchen(p, beta = beta, lambda = lambda,
                lower.tail = lower.tail, log.p = log.p)
}

#' @rdname chen
#' @export
rchen <- function(n, beta, lambda) {
  .chk_pos(beta, lambda)
  reliaR::qchen(stats::runif(n), beta = beta, lambda = lambda)
}

#' @rdname chen
#' @export
schen <- function(x, beta, lambda) reliaR::schen(x, beta = beta, lambda = lambda)

#' @rdname chen
#' @export
hchen <- function(x, beta, lambda) reliaR::hchen(x, beta = beta, lambda = lambda)

#' @rdname chen
#' @export
hra.chen <- function(x, beta, lambda) {
  reliaR::hra.chen(x, beta = beta, lambda = lambda)
}

#' @rdname chen
#' @export
crf.chen <- function(x, t = 0, beta, lambda) {
  reliaR::crf.chen(x, t = t, beta = beta, lambda = lambda)
}
