#' Bivariate Chen distribution: density, survival, hazard and random
#' generation
#'
#' Probability density, survival function, hazard function and random
#' generation for the Marshall-Olkin type Bivariate Chen distribution
#' \eqn{BvCh(\alpha_1,\alpha_2,\alpha_3,\beta)}.
#'
#' @details Let \eqn{U_i \sim Chen(\alpha_i,\beta)}, \eqn{i=1,2,3}, be
#'   independent and set \eqn{Z_1=\min(U_1,U_3)}, \eqn{Z_2=\min(U_2,U_3)}.
#'   Then the joint survival function is
#'   \deqn{S(z_1,z_2) = \exp\left\{\alpha_1\left(1-e^{z_1^{\beta}}\right)
#'     \right\}\exp\left\{\alpha_2\left(1-e^{z_2^{\beta}}\right)\right\}
#'     \exp\left\{\alpha_3\left(1-e^{w^{\beta}}\right)\right\},}
#'   where \eqn{w=\max(z_1,z_2)}. Each factor is obtained from
#'   \code{\link[reliaR]{schen}}. The joint density is piecewise,
#'   \deqn{f(z_1,z_2) = \begin{cases} f_1(z_1,z_2) & z_1<z_2\\
#'     f_2(z_1,z_2) & z_1>z_2\\ f_0(z) & z_1=z_2=z,\end{cases}}
#'   with
#'   \deqn{f_1=\alpha_1(\alpha_2+\alpha_3)\beta^2 z_1^{\beta-1}z_2^{\beta-1}
#'     \exp\left\{z_1^{\beta}+z_2^{\beta}+\alpha_1(1-e^{z_1^{\beta}})
#'     +(\alpha_2+\alpha_3)(1-e^{z_2^{\beta}})\right\},}
#'   \deqn{f_2=\alpha_2(\alpha_1+\alpha_3)\beta^2 z_1^{\beta-1}z_2^{\beta-1}
#'     \exp\left\{z_1^{\beta}+z_2^{\beta}+\alpha_2(1-e^{z_2^{\beta}})
#'     +(\alpha_1+\alpha_3)(1-e^{z_1^{\beta}})\right\},}
#'   \deqn{f_0(z) = \alpha_3\beta z^{\beta-1}
#'     \exp\left\{z^{\beta}+(\alpha_1+\alpha_2+\alpha_3)(1-e^{z^{\beta}})\right\}.}
#'
#'   The mass carried by the three regions is respectively
#'   \eqn{\alpha_1/(\alpha_1+\alpha_2+\alpha_3)},
#'   \eqn{\alpha_2/(\alpha_1+\alpha_2+\alpha_3)} and
#'   \eqn{\alpha_3/(\alpha_1+\alpha_2+\alpha_3)}, so the density integrates to
#'   one. The hazard function is the ratio \eqn{f/S}, that is
#'   \deqn{h(z_1,z_2)=\begin{cases}
#'     \alpha_1(\alpha_2+\alpha_3)\beta^2 z_1^{\beta-1}z_2^{\beta-1}
#'       e^{z_1^{\beta}+z_2^{\beta}} & z_1<z_2\\
#'     \alpha_2(\alpha_1+\alpha_3)\beta^2 z_1^{\beta-1}z_2^{\beta-1}
#'       e^{z_1^{\beta}+z_2^{\beta}} & z_1>z_2\\
#'     \alpha_3\beta z^{\beta-1}e^{z^{\beta}} & z_1=z_2=z.\end{cases}}
#'
#'   The model is a mixture (Proposition 1) of an absolutely continuous part
#'   of weight \eqn{(\alpha_1+\alpha_2)/(\alpha_1+\alpha_2+\alpha_3)} and a
#'   singular part of weight \eqn{\alpha_3/(\alpha_1+\alpha_2+\alpha_3)}
#'   supported on the diagonal; see \code{\link{dbvch_ac}} and
#'   \code{\link{dbvch_sing}}.
#'
#' @param z1,z2 numeric vectors of first and second component.
#' @param alpha1,alpha2,alpha3,beta positive parameters.
#' @param n number of observations.
#' @param log logical; if \code{TRUE} the log density is returned.
#'
#' @return A numeric vector, or for \code{rbvch} a two column matrix with
#'   columns \code{z1} and \code{z2}.
#'
#' @seealso \code{\link{sbvch}}, \code{\link{hbvch}}, \code{\link{rbvch}},
#'   \code{\link{fitbvch}}, \code{\link{plot_surface_bvch}}
#'
#' @examples
#' dbvch(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#' sbvch(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#' hbvch(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#' set.seed(1)
#' rbvch(5, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#'
#' @name bvch
NULL
#' @rdname bvch
#' @export
dbvch <- function(z1, z2, alpha1, alpha2, alpha3, beta, log = FALSE) {
  .chk_pos(alpha1, alpha2, alpha3, beta)
  r <- .recycle(z1, z2, alpha1, alpha2, alpha3, beta)
  out <- .dbvch_log(r[[1]], r[[2]], r[[3]], r[[4]], r[[5]], r[[6]])
  if (log) out else exp(out)
}

#' @keywords internal
#' @noRd
.dbvch_log <- function(z1, z2, alpha1, alpha2, alpha3, beta) {
  neg <- (z1 < 0) | (z2 < 0)
  zp <- pmax(z1, 0)
  zp2 <- pmax(z2, 0)
  y1 <- zp^beta
  y2 <- zp2^beta
  yw <- pmax(y1, y2)
  lt <- alpha1 * (1 - exp(y1)) + alpha2 * (1 - exp(y2)) +
    alpha3 * (1 - exp(yw))
  base <- 2 * log(beta) + y1 + y2 + lt
  lz <- rep(0, length(zp))
  ne <- beta != 1
  lz[ne] <- (beta[ne] - 1) * (log(zp[ne]) + log(zp2[ne]))
  below <- zp < zp2
  above <- zp > zp2
  eq <- zp == zp2
  out <- numeric(length(zp))
  out[below] <- log(alpha1[below]) + log(alpha2[below] + alpha3[below]) +
    base[below] + lz[below]
  out[above] <- log(alpha2[above]) + log(alpha1[above] + alpha3[above]) +
    base[above] + lz[above]
  if (any(eq)) {
    le <- (beta[eq] - 1) * log(zp[eq])
    le[beta[eq] == 1] <- 0
    out[eq] <- log(alpha3[eq]) + log(beta[eq]) + le + y1[eq] +
      (alpha1[eq] + alpha2[eq] + alpha3[eq]) * (1 - exp(y1[eq]))
  }
  out[neg] <- -Inf
  out
}

#' @rdname bvch
#' @export
sbvch <- function(z1, z2, alpha1, alpha2, alpha3, beta) {
  .chk_pos(alpha1, alpha2, alpha3, beta)
  r <- .recycle(z1, z2, alpha1, alpha2, alpha3, beta)
  w <- pmax(pmax(r[[1]], 0), pmax(r[[2]], 0))
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]
  exp(a1 * (1 - exp(pmax(r[[1]], 0)^b)) +
    a2 * (1 - exp(pmax(r[[2]], 0)^b)) +
    a3 * (1 - exp(w^b)))
}

#' @rdname bvch
#' @export
hbvch <- function(z1, z2, alpha1, alpha2, alpha3, beta) {
  .chk_pos(alpha1, alpha2, alpha3, beta)
  r <- .recycle(z1, z2, alpha1, alpha2, alpha3, beta)
  zp <- pmax(r[[1]], 0); zp2 <- pmax(r[[2]], 0)
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]
  y1 <- zp^b; y2 <- zp2^b
  lz <- rep(0, length(zp))
  ne <- b != 1
  lz[ne] <- (b[ne] - 1) * (log(zp[ne]) + log(zp2[ne]))
  base <- 2 * log(b) + y1 + y2 + lz
  below <- zp < zp2
  above <- zp > zp2
  eq <- zp == zp2
  out <- numeric(length(zp))
  out[below] <- exp(log(a1[below]) + log(a2[below] + a3[below]) + base[below])
  out[above] <- exp(log(a2[above]) + log(a1[above] + a3[above]) + base[above])
  if (any(eq)) {
    le <- (b[eq] - 1) * log(zp[eq])
    le[b[eq] == 1] <- 0
    out[eq] <- exp(log(a3[eq]) + log(b[eq]) + le + y1[eq])
  }
  out[zp < 0 | zp2 < 0] <- 0
  out
}

#' @rdname bvch
#' @export
rbvch <- function(n, alpha1, alpha2, alpha3, beta) {
  .chk_pos(alpha1, alpha2, alpha3, beta)
  n <- as.integer(n)
  u1 <- .rChen(n, beta, alpha1)
  u2 <- .rChen(n, beta, alpha2)
  u3 <- .rChen(n, beta, alpha3)
  cbind(z1 = pmin(u1, u3), z2 = pmin(u2, u3))
}

#' Absolutely continuous and singular components of the BvCh density
#'
#' Decomposition of the Bivariate Chen density into the two components of the
#' mixture representation
#' \eqn{f = \frac{\alpha_1+\alpha_2}{\alpha_1+\alpha_2+\alpha_3}f_a +
#' \frac{\alpha_3}{\alpha_1+\alpha_2+\alpha_3}f_s}.
#'
#' @details \code{dbvch_ac} is the properly normalised absolutely continuous
#'   component of Proposition 1; it integrates to one over the plane, and
#'   \code{dbvch_sing} integrates to one along the diagonal. Consequently
#'   \code{dbvch} is recovered exactly as
#'   \code{(alpha1 + alpha2)/(alpha1 + alpha2 + alpha3) * dbvch_ac(...)} plus
#'   \code{alpha3/(alpha1 + alpha2 + alpha3) * dbvch_sing(...)}.
#'
#' @inheritParams bvch
#' @return A numeric vector.
#' @seealso \code{\link{dbvch}}
#' @examples
#' a <- 1.5; b2 <- 2.5; c3 <- 0.5; bt <- 1.8
#' zz1 <- c(0.4, 0.9, 0.7); zz2 <- c(0.6, 0.4, 0.7)
#' ac <- dbvch_ac(zz1, zz2, a, b2, c3, bt)
#' sg <- dbvch_sing(zz1, zz2, a, b2, c3, bt)
#' k <- (a + b2 + c3)
#' mix <- (a + b2) / k * ac + c3 / k * sg
#' all.equal(mix, dbvch(zz1, zz2, a, b2, c3, bt))
#'
#' @name bvch_components
NULL

#' @rdname bvch_components
#' @export
dbvch_ac <- function(z1, z2, alpha1, alpha2, alpha3, beta) {
  .chk_pos(alpha1, alpha2, alpha3, beta)
  r <- .recycle(z1, z2, alpha1, alpha2, alpha3, beta)
  zp <- pmax(r[[1]], 0); zp2 <- pmax(r[[2]], 0)
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]
  k <- (a1 + a2 + a3) / (a1 + a2)
  f <- dbvch(zp, zp2, a1, a2, a3, b)
  f[zp == zp2] <- 0
  k * f
}

#' @rdname bvch_components
#' @export
dbvch_sing <- function(z1, z2, alpha1, alpha2, alpha3, beta) {
  .chk_pos(alpha1, alpha2, alpha3, beta)
  r <- .recycle(z1, z2, alpha1, alpha2, alpha3, beta)
  zp <- pmax(r[[1]], 0); zp2 <- pmax(r[[2]], 0)
  a1 <- r[[3]]; a2 <- r[[4]]; a3 <- r[[5]]; b <- r[[6]]
  out <- reliaR::dchen(pmax(zp, .Machine$double.xmin), b, a1 + a2 + a3)
  out[zp != zp2] <- 0
  out[(r[[1]] < 0) | (r[[2]] < 0)] <- 0
  out
}
