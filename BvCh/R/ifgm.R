#' Inverse FGM copula
#'
#' The Inverse FGM (IFGM) copula, its density, dependence measures and a
#' random generator.
#'
#' @details For \eqn{0\le\gamma,\omega\le 1} the Inverse FGM copula is
#' \deqn{C(u,v)=\frac{uv}{\sqrt{1+\gamma(1-u)(1-v)+\omega uv(1-u)(1-v)}}.}
#' Writing \eqn{D=1+\gamma(1-u)(1-v)+\omega uv(1-u)(1-v)} and
#' \eqn{P=uv}, its density has the closed form
#' \deqn{c(u,v)=D^{-1/2}
#'   -\frac{v D_v + u D_u}{2D^{3/2}}
#'   +\frac{3P D_u D_v}{4D^{5/2}}
#'   -\frac{P D_{uv}}{2D^{3/2}},}
#' with \eqn{D_u=-\gamma(1-v)+\omega(1-2u)v(1-v)},
#' \eqn{D_v=-\gamma(1-u)+\omega u(1-2v)(1-u)} and
#' \eqn{D_{uv}=\gamma+\omega(1-2u)(1-2v)}.
#'
#' Kendall's tau and Spearman's rho are computed as
#' \deqn{\tau=4\int_0^1\!\int_0^1 C(u,v)c(u,v)\,du\,dv-1,\qquad
#'   \rho=12\int_0^1\!\int_0^1 uv\,c(u,v)\,du\,dv-3.}
#' Both are obtained by numerical quadrature and reproduce Monte Carlo
#' estimates from \code{rifgm}.
#'
#' @section Sign of the dependence:
#' With the \eqn{+\,\gamma(1-u)(1-v)+\omega uv(1-u)(1-v)} denominator given
#' above, the copula lies \emph{below} the independence copula
#' \eqn{uv} (\eqn{C(0.5,0.5)=0.2306<0.25} at \eqn{\gamma=0.6,\omega=0.4}),
#' so \eqn{\tau} and \eqn{\rho} are non-positive for \eqn{\gamma,\omega\ge0}:
#' numerically \eqn{\tau=-0.0699} and \eqn{\rho=-0.1049} at
#' \eqn{\gamma=0.6,\omega=0.4}. This is the usual behaviour of an
#' \emph{inverse} FGM type copula. Both measures are exactly zero at
#' \eqn{\gamma=\omega=0}. Note that the summary document this package was
#' written from quotes \eqn{\tau\approx 0.156>0} for the same parameters,
#' which is inconsistent with the copula formula it also quotes; the formula
#' has been implemented as given, and \code{tau_ifgm}/\code{rho_ifgm} are
#' provided so the dependence can be checked directly.
#'
#' @param u,v numeric vectors of copula arguments in \eqn{[0,1]}.
#' @param gamma,omega copula parameters, each in \eqn{[0,1]}.
#' @param n number of observations.
#'
#' @return \code{cifgm} gives the copula, \code{cifgmdens} its density,
#'   \code{rifgm} a two column matrix of copula draws, and \code{tau_ifgm}
#'   and \code{rho_ifgm} scalar dependence measures.
#'
#' @seealso \code{\link{difgmchen}}, \code{\link{hifgmchen}}
#'
#' @examples
#' cifgm(0.3, 0.7, gamma = 0.6, omega = 0.4)
#' cifgmdens(0.3, 0.7, gamma = 0.6, omega = 0.4)
#' tau_ifgm(0.6, 0.4)
#' rho_ifgm(0.6, 0.4)
#' set.seed(1)
#' head(rifgm(5, gamma = 0.6, omega = 0.4))
#'
#' @name ifgm
NULL

#' @keywords internal
#' @noRd
.ifgm_D <- function(u, v, gamma, omega) {
  1 + gamma * (1 - u) * (1 - v) + omega * u * v * (1 - u) * (1 - v)
}

#' @rdname ifgm
#' @export
cifgm <- function(u, v, gamma, omega) {
  if (any(gamma < 0 | gamma > 1) || any(omega < 0 | omega > 1))
    stop("gamma and omega must lie in [0, 1].", call. = FALSE)
  r <- .recycle(u, v, gamma, omega)
  u <- pmin(pmax(r[[1]], 0), 1); v <- pmin(pmax(r[[2]], 0), 1)
  r[[1]] * r[[2]] / sqrt(.ifgm_D(u, v, r[[3]], r[[4]]))
}

#' @rdname ifgm
#' @export
cifgmdens <- function(u, v, gamma, omega) {
  if (any(gamma < 0 | gamma > 1) || any(omega < 0 | omega > 1))
    stop("gamma and omega must lie in [0, 1].", call. = FALSE)
  r <- .recycle(u, v, gamma, omega)
  u <- pmin(pmax(r[[1]], 0), 1); v <- pmin(pmax(r[[2]], 0), 1)
  gg <- r[[3]]; om <- r[[4]]
  D <- .ifgm_D(u, v, gg, om)
  P <- u * v
  Du <- -gg * (1 - v) + om * (1 - 2 * u) * v * (1 - v)
  Dv <- -gg * (1 - u) + om * u * (1 - 2 * v) * (1 - u)
  Duv <- gg + om * (1 - 2 * u) * (1 - 2 * v)
  D^(-0.5) - (v * Dv + u * Du) / (2 * D^1.5) +
    3 * P * Du * Dv / (4 * D^2.5) - P * Duv / (2 * D^1.5)
}

#' Conditional distribution of the second argument given the first
#' @keywords internal
#' @noRd
.cifgm_cond <- function(v, u, gamma, omega) {
  ## dC/du evaluated at (u, v); C(u, 1) = u, so this is the conditional cdf
  D <- .ifgm_D(u, v, gamma, omega)
  Du <- -gamma * (1 - v) + omega * (1 - 2 * u) * v * (1 - v)
  v * D^(-0.5) - (u * v / 2) * D^(-1.5) * Du
}

#' @rdname ifgm
#' @export
rifgm <- function(n, gamma, omega) {
  if (any(gamma < 0 | gamma > 1) || any(omega < 0 | omega > 1))
    stop("gamma and omega must lie in [0, 1].", call. = FALSE)
  n <- as.integer(n)
  u <- stats::runif(n)
  w <- stats::runif(n)
  ## v is obtained by inverting the conditional cdf .cifgm_cond, which is
  ## increasing from 0 to 1 on [0, 1].  A vectorised bisection keeps this
  ## exact to machine precision and linear in n.
  lo <- rep(0, n)
  hi <- rep(1, n)
  for (i in seq_len(64L)) {
    mid <- (lo + hi) / 2
    below <- .cifgm_cond(mid, u, gamma, omega) < w
    lo[below] <- mid[below]
    hi[!below] <- mid[!below]
  }
  v <- (lo + hi) / 2
  cbind(u = u, v = pmin(pmax(v, 0), 1))
}

#' @rdname ifgm
#' @export
tau_ifgm <- function(gamma, omega) {
  f <- function(u) {
    vapply(u, function(uu) {
      stats::integrate(function(vv)
        cifgm(uu, vv, gamma, omega) * cifgmdens(uu, vv, gamma, omega),
        0, 1, rel.tol = 1e-10)$value
    }, numeric(1))
  }
  4 * stats::integrate(f, 0, 1, rel.tol = 1e-10)$value - 1
}

#' @rdname ifgm
#' @export
rho_ifgm <- function(gamma, omega) {
  f <- function(u) {
    vapply(u, function(uu) {
      stats::integrate(function(vv)
        uu * vv * cifgmdens(uu, vv, gamma, omega), 0, 1,
        rel.tol = 1e-10)$value
    }, numeric(1))
  }
  12 * stats::integrate(f, 0, 1, rel.tol = 1e-10)$value - 3
}
#' Bivariate IFGM-Chen distribution
#'
#' Reliability function, density, hazard function and random generation for
#' the Bivariate Inverse FGM Chen model, obtained by combining the Inverse
#' FGM copula with Chen marginals.
#'
#' @details The marginals are Chen distributions
#' \eqn{F_X(x)=1-e^{-\theta_1(e^{x^{\beta_1}}-1)}} and
#' \eqn{F_Y(y)=1-e^{-\theta_2(e^{y^{\beta_2}}-1)}}, which are exactly
#' \code{reliaR::pchen(x, beta1, theta1)} and
#' \code{reliaR::pchen(y, beta2, theta2)}. With \eqn{C} the Inverse FGM
#' copula of \code{\link{cifgm}} and \eqn{c} its density, the joint
#' reliability function and density are
#' \deqn{R(x,y)=C(F_X(x),F_Y(y))-F_X(x)-F_Y(y)+1,}
#' \deqn{f(x,y)=c(F_X(x),F_Y(y))\,f_X(x)\,f_Y(y),}
#' where \eqn{f_X = dchen(x,\beta_1,\theta_1)} and
#' \eqn{f_Y = dchen(y,\beta_2,\theta_2)}, where \code{dchen} is
#' \code{\link[reliaR]{dchen}}. The hazard function is
#' the ratio \eqn{f/R}.
#'
#' @section Remarks:
#' Unlike the Marshall-Olkin type models, this construction is absolutely
#' continuous: it has no mass on the diagonal and therefore no singular
#' component. The survival function is
#' \code{1 - pchen(x, beta1, theta1) - pchen(y, beta2, theta2) + R(x, y)}
#' and can be negative by cancellation in the far tail.
#'
#' @param x,y numeric vectors of first and second component.
#' @param theta1,beta1,theta2,beta2 positive marginal parameters.
#' @param gamma,omega copula parameters in \eqn{[0,1]}.
#' @param n number of observations.
#' @param log logical; if \code{TRUE} the log density is returned.
#'
#' @return A numeric vector, or for \code{rifgmchen} a two column matrix
#'   with columns \code{x} and \code{y}.
#'
#' @seealso \code{\link{cifgm}}, \code{\link{fitifgmchen}},
#'   \code{\link{plot_reliability_ifgmchen}}
#'
#' @examples
#' difgmchen(0.5, 0.9, theta1 = 1.2, beta1 = 1.3, theta2 = 1.5,
#'           beta2 = 1.1, gamma = 0.6, omega = 0.4)
#' relifgmchen(0.5, 0.9, theta1 = 1.2, beta1 = 1.3, theta2 = 1.5,
#'             beta2 = 1.1, gamma = 0.6, omega = 0.4)
#' hifgmchen(0.5, 0.9, theta1 = 1.2, beta1 = 1.3, theta2 = 1.5,
#'           beta2 = 1.1, gamma = 0.6, omega = 0.4)
#' set.seed(1)
#' head(rifgmchen(5, theta1 = 1.2, beta1 = 1.3, theta2 = 1.5,
#'                beta2 = 1.1, gamma = 0.6, omega = 0.4))
#'
#' @name ifgmchen
NULL

#' @rdname ifgmchen
#' @export
relifgmchen <- function(x, y, theta1, beta1, theta2, beta2, gamma, omega) {
  .chk_pos(theta1, beta1, theta2, beta2)
  r <- .recycle(x, y, theta1, beta1, theta2, beta2, gamma, omega)
  Fx <- .pChen(pmax(r[[1]], 0), r[[4]], r[[3]])
  Fy <- .pChen(pmax(r[[2]], 0), r[[6]], r[[5]])
  cifgm(Fx, Fy, r[[7]], r[[8]]) - Fx - Fy + 1
}

#' @rdname ifgmchen
#' @export
difgmchen <- function(x, y, theta1, beta1, theta2, beta2, gamma, omega,
                      log = FALSE) {
  .chk_pos(theta1, beta1, theta2, beta2)
  r <- .recycle(x, y, theta1, beta1, theta2, beta2, gamma, omega)
  xp <- pmax(r[[1]], 0); yp <- pmax(r[[2]], 0)
  t1 <- r[[3]]; b1 <- r[[4]]; t2 <- r[[5]]; b2 <- r[[6]]
  gg <- r[[7]]; om <- r[[8]]
  Fx <- .pChen(xp, b1, t1)
  Fy <- .pChen(yp, b2, t2)
  out <- cifgmdens(Fx, Fy, gg, om) * reliaR::dchen(pmax(xp, .Machine$double.xmin), b1, t1) *
    reliaR::dchen(pmax(yp, .Machine$double.xmin), b2, t2)
  out[(r[[1]] < 0) | (r[[2]] < 0)] <- 0
  if (log) log(out) else out
}

#' @rdname ifgmchen
#' @export
hifgmchen <- function(x, y, theta1, beta1, theta2, beta2, gamma, omega) {
  .chk_pos(theta1, beta1, theta2, beta2)
  R <- relifgmchen(x, y, theta1, beta1, theta2, beta2, gamma, omega)
  f <- difgmchen(x, y, theta1, beta1, theta2, beta2, gamma, omega)
  out <- f / R
  out[!is.finite(out)] <- 0
  out
}

#' @rdname ifgmchen
#' @export
rifgmchen <- function(n, theta1, beta1, theta2, beta2, gamma, omega) {
  .chk_pos(theta1, beta1, theta2, beta2)
  n <- as.integer(n)
  uv <- rifgm(n, gamma, omega)
  cbind(x = reliaR::qchen(uv[, 1], beta1, theta1),
        y = reliaR::qchen(uv[, 2], beta2, theta2))
}
