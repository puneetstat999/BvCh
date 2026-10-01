#' Bivariate Extended Chen distribution: density, survival, hazard and
#' random generation
#'
#' Probability density, survival function, hazard function and random
#' generation for the Marshall-Olkin type Bivariate Extended Chen
#' distribution \eqn{BEC(\alpha_1,\alpha_2,\alpha_3,\beta,\lambda)}.
#'
#' @details With \eqn{U_i \sim EC(\alpha_i,\beta,\lambda)} independent and
#'   \eqn{X_1=\min(U_1,U_3)}, \eqn{X_2=\min(U_2,U_3)}, the joint survival
#'   function is
#'   \deqn{S(x_1,x_2)=\begin{cases}
#'     S_{EC}(x_1;\alpha_1,\beta,\lambda)\,
#'       S_{EC}(x_2;\alpha_2+\alpha_3,\beta,\lambda) & x_2>x_1\\
#'     S_{EC}(x_1;\alpha_1+\alpha_3,\beta,\lambda)\,
#'       S_{EC}(x_2;\alpha_2,\beta,\lambda) & x_1>x_2\\
#'     S_{EC}(x;\alpha_1+\alpha_2+\alpha_3,\beta,\lambda) & x_1=x_2=x,\end{cases}}
#'   where \eqn{S_{EC}(t;\alpha,\beta,\lambda)=[1+\lambda(e^{t^\beta}-1)]^{-\alpha}}
#'   (see \code{\link{sech}}). The joint density is
#'   \deqn{f(x_1,x_2)=\begin{cases} f_1 & x_2>x_1\\ f_2 & x_1>x_2\\
#'     f_3 & x_1=x_2=x,\end{cases}}
#'   with
#'   \deqn{f_1=\alpha_1(\alpha_2+\alpha_3)\lambda^2\beta^2
#'     x_1^{\beta-1}x_2^{\beta-1}e^{x_1^{\beta}+x_2^{\beta}}
#'     [1+\lambda(e^{x_1^{\beta}}-1)]^{-(\alpha_1+1)}
#'     [1+\lambda(e^{x_2^{\beta}}-1)]^{-(\alpha_2+\alpha_3+1)},}
#'   \deqn{f_2=\alpha_2(\alpha_1+\alpha_3)\lambda^2\beta^2
#'     x_1^{\beta-1}x_2^{\beta-1}e^{x_1^{\beta}+x_2^{\beta}}
#'     [1+\lambda(e^{x_1^{\beta}}-1)]^{-(\alpha_1+\alpha_3+1)}
#'     [1+\lambda(e^{x_2^{\beta}}-1)]^{-(\alpha_2+1)},}
#'   \deqn{f_3(x)=\alpha_3\beta\lambda x^{\beta-1}e^{x^{\beta}}
#'     [1+\lambda(e^{x^{\beta}}-1)]^{-(\alpha_1+\alpha_2+\alpha_3+1)}.}
#'   The three regions carry mass \eqn{\alpha_i/(\alpha_1+\alpha_2+\alpha_3)},
#'   so the density integrates to one. The hazard function is \eqn{f/S},
#'   \deqn{h(x_1,x_2)=\begin{cases}
#'     \dfrac{\alpha_1(\alpha_2+\alpha_3)\lambda^2\beta^2x_1^{\beta-1}
#'       x_2^{\beta-1}e^{x_1^{\beta}+x_2^{\beta}}}
#'       {[1+\lambda(e^{x_1^{\beta}}-1)][1+\lambda(e^{x_2^{\beta}}-1)]}
#'       & x_2>x_1\\[2ex]
#'     \dfrac{\alpha_2(\alpha_1+\alpha_3)\lambda^2\beta^2x_1^{\beta-1}
#'       x_2^{\beta-1}e^{x_1^{\beta}+x_2^{\beta}}}
#'       {[1+\lambda(e^{x_1^{\beta}}-1)][1+\lambda(e^{x_2^{\beta}}-1)]}
#'       & x_1>x_2\\[2ex]
#'     \dfrac{\alpha_3\beta\lambda x^{\beta-1}e^{x^{\beta}}}
#'       {1+\lambda(e^{x^{\beta}}-1)} & x_1=x_2=x.\end{cases}}
#'   The componentwise hazard gradients are available separately through
#'   \code{\link{hbec1}} and \code{\link{hbec2}}.
#'
#' @param x1,x2 numeric vectors of first and second component.
#' @param alpha1,alpha2,alpha3,beta,lambda positive parameters.
#' @param n number of observations.
#' @param log logical; if \code{TRUE} the log density is returned.
#'
#' @return A numeric vector, or for \code{rbec} a two column matrix with
#'   columns \code{x1} and \code{x2}.
#'
#' @seealso \code{\link{sbec}}, \code{\link{hbec}}, \code{\link{rbec}},
#'   \code{\link{fitbec}}, \code{\link{sech}}
#'
#' @examples
#' dbec(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8, lambda = 2)
#' sbec(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8, lambda = 2)
#' hbec(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8, lambda = 2)
#' set.seed(1)
#' rbec(5, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8, lambda = 2)
#'
#' @name bec
NULL

#' @keywords internal
#' @noRd
.lech <- function(x, alpha, beta, lambda) {
  g <- .log1p_lam_expm1(x^beta, lambda)
  -alpha * g
}

#' @keywords internal
#' @noRd
.hech <- function(x, alpha, beta, lambda) {
  x <- pmax(x, 0)
  y <- x^beta
  lz <- rep(0, length(x))
  ne <- beta != 1
  lz[ne] <- (beta[ne] - 1) * log(x[ne])
  out <- exp(log(alpha) + log(lambda) + log(beta) + lz + y -
    .log1p_lam_expm1(y, lambda))
  out
}

#' @rdname bec
#' @export
dbec <- function(x1, x2, alpha1, alpha2, alpha3, beta, lambda, log = FALSE) {
  .chk_pos(alpha1, alpha2, alpha3, beta, lambda)
  r <- .recycle(x1, x2, alpha1, alpha2, alpha3, beta, lambda)
  xp <- pmax(r[[1]], 0); xp2 <- pmax(r[[2]], 0)
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]; lam <- r[[7]]
  y1 <- xp^b; y2 <- xp2^b
  g1 <- .log1p_lam_expm1(y1, lam)
  g2 <- .log1p_lam_expm1(y2, lam)
  lz <- rep(0, length(xp))
  ne <- b != 1
  lz[ne] <- (b[ne] - 1) * (log(xp[ne]) + log(xp2[ne]))
  base <- 2 * log(b) + 2 * log(lam) + lz + y1 + y2
  below <- xp < xp2
  above <- xp > xp2
  eq <- xp == xp2
  out <- numeric(length(xp))
  out[below] <- log(a1[below]) + log(a2[below] + a3[below]) + base[below] -
    (a1[below] + 1) * g1[below] - (a2[below] + a3[below] + 1) * g2[below]
  out[above] <- log(a2[above]) + log(a1[above] + a3[above]) + base[above] -
    (a1[above] + a3[above] + 1) * g1[above] - (a2[above] + 1) * g2[above]
  if (any(eq)) {
    le <- (b[eq] - 1) * log(xp[eq])
    le[b[eq] == 1] <- 0
    out[eq] <- log(a3[eq]) + log(b[eq]) + log(lam[eq]) + le + y1[eq] -
      (a1[eq] + a2[eq] + a3[eq] + 1) * g1[eq]
  }
  out[(r[[1]] < 0) | (r[[2]] < 0)] <- -Inf
  if (log) out else exp(out)
}

#' @rdname bec
#' @export
sbec <- function(x1, x2, alpha1, alpha2, alpha3, beta, lambda) {
  .chk_pos(alpha1, alpha2, alpha3, beta, lambda)
  r <- .recycle(x1, x2, alpha1, alpha2, alpha3, beta, lambda)
  xp <- pmax(r[[1]], 0); xp2 <- pmax(r[[2]], 0)
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]; lam <- r[[7]]
  below <- xp < xp2
  above <- xp > xp2
  eq <- xp == xp2
  g1 <- .log1p_lam_expm1(xp^b, lam)
  g2 <- .log1p_lam_expm1(xp2^b, lam)
  out <- numeric(length(xp))
  out[below] <- -a1[below] * g1[below] - (a2[below] + a3[below]) * g2[below]
  out[above] <- -(a1[above] + a3[above]) * g1[above] - a2[above] * g2[above]
  out[eq] <- -(a1[eq] + a2[eq] + a3[eq]) * g1[eq]
  exp(out)
}

#' @rdname bec
#' @export
hbec <- function(x1, x2, alpha1, alpha2, alpha3, beta, lambda) {
  .chk_pos(alpha1, alpha2, alpha3, beta, lambda)
  r <- .recycle(x1, x2, alpha1, alpha2, alpha3, beta, lambda)
  xp <- pmax(r[[1]], 0); xp2 <- pmax(r[[2]], 0)
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]; lam <- r[[7]]
  y1 <- xp^b; y2 <- xp2^b
  g1 <- .log1p_lam_expm1(y1, lam)
  g2 <- .log1p_lam_expm1(y2, lam)
  lz <- rep(0, length(xp))
  ne <- b != 1
  lz[ne] <- (b[ne] - 1) * (log(xp[ne]) + log(xp2[ne]))
  base <- 2 * log(b) + 2 * log(lam) + lz + y1 + y2 - g1 - g2
  below <- xp < xp2
  above <- xp > xp2
  eq <- xp == xp2
  out <- numeric(length(xp))
  out[below] <- exp(log(a1[below]) + log(a2[below] + a3[below]) + base[below])
  out[above] <- exp(log(a2[above]) + log(a1[above] + a3[above]) + base[above])
  if (any(eq)) {
    le <- (b[eq] - 1) * log(xp[eq])
    le[b[eq] == 1] <- 0
    out[eq] <- exp(log(a3[eq]) + log(b[eq]) + log(lam[eq]) + le + y1[eq] - g1[eq])
  }
  out[(r[[1]] < 0) | (r[[2]] < 0)] <- 0
  out
}

#' @rdname bec
#' @export
rbec <- function(n, alpha1, alpha2, alpha3, beta, lambda) {
  .chk_pos(alpha1, alpha2, alpha3, beta, lambda)
  n <- as.integer(n)
  u1 <- rech(n, alpha1, beta, lambda)
  u2 <- rech(n, alpha2, beta, lambda)
  u3 <- rech(n, alpha3, beta, lambda)
  cbind(x1 = pmin(u1, u3), x2 = pmin(u2, u3))
}

#' Componentwise hazard gradients of the BEC model
#'
#' The hazard gradients of the Bivariate Extended Chen model, that is the
#' hazard of the \eqn{i}-th marginal conditional on the other component.
#'
#' @details
#' \deqn{h_{X_1}(x_1,x_2)=\begin{cases}
#'   \dfrac{\alpha_1\beta\lambda x_1^{\beta-1}e^{x_1^{\beta}}}
#'     {1+\lambda(e^{x_1^{\beta}}-1)} & x_2>x_1\\
#'   \dfrac{(\alpha_1+\alpha_3)\beta\lambda x_1^{\beta-1}e^{x_1^{\beta}}}
#'     {1+\lambda(e^{x_1^{\beta}}-1)} & x_1>x_2\\
#'   \dfrac{(\alpha_1+\alpha_2+\alpha_3)\beta\lambda x_1^{\beta-1}
#'     e^{x_1^{\beta}}}{1+\lambda(e^{x_1^{\beta}}-1)} & x_1=x_2,\end{cases}}
#' with \eqn{h_{X_2}} obtained by interchanging the roles of \eqn{x_1},
#' \eqn{x_2} and \eqn{\alpha_1}, \eqn{\alpha_2}.
#'
#' @inheritParams bec
#' @return A numeric vector.
#' @seealso \code{\link{hbec}}
#' @examples
#' hbec1(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8, lambda = 2)
#' hbec2(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8, lambda = 2)
#'
#' @name hbec_grad
NULL

#' @rdname hbec_grad
#' @export
hbec1 <- function(x1, x2, alpha1, alpha2, alpha3, beta, lambda) {
  .chk_pos(alpha1, alpha2, alpha3, beta, lambda)
  r <- .recycle(x1, x2, alpha1, alpha2, alpha3, beta, lambda)
  xp <- pmax(r[[1]], 0); xp2 <- pmax(r[[2]], 0)
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]; lam <- r[[7]]
  below <- xp < xp2
  above <- xp > xp2
  eq <- xp == xp2
  al <- ifelse(below, a1, ifelse(above, a1 + a3, a1 + a2 + a3))
  h <- .hech(xp, al, b, lam)
  h[(r[[1]] < 0) | (r[[2]] < 0)] <- 0
  h
}

#' @rdname hbec_grad
#' @export
hbec2 <- function(x1, x2, alpha1, alpha2, alpha3, beta, lambda) {
  .chk_pos(alpha1, alpha2, alpha3, beta, lambda)
  r <- .recycle(x1, x2, alpha1, alpha2, alpha3, beta, lambda)
  xp <- pmax(r[[1]], 0); xp2 <- pmax(r[[2]], 0)
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]; lam <- r[[7]]
  below <- xp < xp2
  above <- xp > xp2
  eq <- xp == xp2
  al <- ifelse(below, a2 + a3, ifelse(above, a2, a1 + a2 + a3))
  h <- .hech(xp2, al, b, lam)
  h[(r[[1]] < 0) | (r[[2]] < 0)] <- 0
  h
}
#' rbec(5, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8, lambda = 2)
#'
#' @name bec
NULL
