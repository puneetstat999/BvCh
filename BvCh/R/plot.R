#' Surface and contour plots for the bivariate Chen models
#'
#' Graphical displays of the density, survival or reliability surfaces of the
#' BvCh, BEC and IFGM-Chen models, using base R graphics only.
#'
#' @details \code{plot_surface_bvch} draws the absolutely continuous part of
#'   the Bivariate Chen density, which is the surface shown for varying
#'   \eqn{\beta} in the source paper: as \eqn{\beta} increases the surface
#'   changes from flat or monotone to peaked near the origin. \code{part}
#'   selects between the absolutely continuous part, the singular part and
#'   the full piecewise density.
#'
#' @param alpha1,alpha2,alpha3,beta,lambda model parameters.
#' @param beta.vec a vector of \eqn{\beta} values, one surface per value.
#' @param theta1,beta1,theta2,beta2 IFGM-Chen marginal parameters.
#' @param gamma,omega IFGM-Chen copula parameters.
#' @param what one of \code{"density"}, \code{"survival"} or
#'   \code{"hazard"}.
#' @param part one of \code{"ac"}, \code{"singular"} or \code{"all"}; only
#'   used when \code{what = "density"} for the BvCh model.
#' @param n.grid resolution of the plotting grid in each direction.
#' @param range maximum value of both axes.
#' @param theta,phi camera angles passed to \code{\link[graphics]{persp}}.
#' @param expand.fac expansion factor for the z axis.
#' @param col.colours colours for the surface: a vector of colour names, a
#'   function such as \code{\link[grDevices]{heat.colors}}, or \code{NULL}
#'   for the default heat palette.
#' @param n.col number of colour bands.
#' @param main optional plot title; a default is built from the parameters.
#' @param lab main title for the contour plot.
#' @param add if \code{TRUE} the contour is drawn on the density surface.
#' @param draw.singular logical; draw the singular support on the diagonal
#'   of a contour plot.
#' @param ... further arguments passed to \code{\link[graphics]{persp}}.
#'
#' @return Invisibly, the matrix of values plotted. \code{plot_surface_bvch}
#'   returns one matrix per element of \code{beta.vec}.
#'
#' @seealso \code{\link{dbvch}}, \code{\link{dbec}},
#'   \code{\link{difgmchen}}
#'
#' @examples
#' plot_surface_bvch(alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5,
#'                   beta.vec = c(0.8, 1.8, 5, 10))
#' plot_surface_bec(alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8,
#'                 lambda = 2)
#' plot_reliability_ifgmchen(theta1 = 1.2, beta1 = 1.3, theta2 = 1.5,
#'                           beta2 = 1.1, gamma = 0.6, omega = 0.4)
#' plot_contour_bvch(alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#'
#' @name plotbvch
NULL

#' @keywords internal
#' @noRd
.grid2 <- function(n, range) seq(0, range, length.out = n)

#' @keywords internal
#' @rdname plotbvch
#' @export
plot_surface_bvch <- function(alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5,
                              beta.vec = c(0.8, 1.8, 5, 10),
                              what = c("density", "survival", "hazard"),
                              part = c("ac", "singular", "all"),
                              n.grid = 61L, range = 3, theta = 30, phi = 25,
                              expand.fac = 0.6, col.colours = NULL, n.col = 100, main = NULL, ...) {
  what <- match.arg(what)
  part <- match.arg(part)
  .chk_pos(alpha1, alpha2, alpha3)
  beta.vec <- beta.vec[beta.vec > 0]
  if (!length(beta.vec))
    stop("'beta.vec' must contain positive values.", call. = FALSE)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  if (length(beta.vec) > 1L)
    graphics::par(mfrow = c(ceiling(length(beta.vec) / 2), 2))
  g <- .grid2(n.grid, range)
  G <- as.matrix(expand.grid(g, g))
  x1 <- G[, 2]; x2 <- G[, 1]
  pal <- .palette(col.colours, n.col)
  out <- vector("list", length(beta.vec))
  for (i in seq_along(beta.vec)) {
    bt <- beta.vec[i]
    z <- switch(what,
      density = switch(part,
        ac = dbvch_ac(x1, x2, alpha1, alpha2, alpha3, bt),
        singular = dbvch_sing(x1, x2, alpha1, alpha2, alpha3, bt),
        all = dbvch(x1, x2, alpha1, alpha2, alpha3, bt)),
      survival = sbvch(x1, x2, alpha1, alpha2, alpha3, bt),
      hazard = hbvch(x1, x2, alpha1, alpha2, alpha3, bt))
    ttl <- if (is.null(main)) {
      paste0("BvCh ", what,
             if (what == "density") paste0(" (", part, ")") else "",
             ": alpha=(", alpha1, ", ", alpha2, ", ", alpha3, "), beta=", bt)
    } else if (length(main) == length(beta.vec)) main[i] else main
    out[[i]] <- .surface(z, n.grid, range, theta, phi, expand.fac, pal, ttl, ...)
  }
  names(out) <- paste0("beta=", beta.vec)
  invisible(out)
}

#' @rdname plotbvch
#' @export
plot_surface_bec <- function(alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5,
                             beta = 1.8, lambda = 2,
                             what = c("density", "survival", "hazard"),
                             n.grid = 61L, range = 3, theta = 30, phi = 25,
                             expand.fac = 0.6, col.colours = NULL, n.col = 100,
                             main = NULL, ...) {
  what <- match.arg(what)
  .chk_pos(alpha1, alpha2, alpha3, beta, lambda)
  g <- .grid2(n.grid, range)
  G <- as.matrix(expand.grid(g, g))
  x1 <- G[, 2]; x2 <- G[, 1]
  z <- switch(what,
    density = dbec(x1, x2, alpha1, alpha2, alpha3, beta, lambda),
    survival = sbec(x1, x2, alpha1, alpha2, alpha3, beta, lambda),
    hazard = hbec(x1, x2, alpha1, alpha2, alpha3, beta, lambda))
  ttl <- if (is.null(main))
    paste0("BEC ", what, ": alpha=(", alpha1, ", ", alpha2, ", ", alpha3,
           "), beta=", beta, ", lambda=", lambda) else main
  .surface(z, n.grid, range, theta, phi, expand.fac,
           .palette(col.colours, n.col), ttl, ...)
}

#' @rdname plotbvch
#' @export
plot_reliability_ifgmchen <- function(theta1 = 1.2, beta1 = 1.3,
                                      theta2 = 1.5, beta2 = 1.1,
                                      gamma = 0.6, omega = 0.4,
                                      n.grid = 61L, range = 3, theta = 30,
                                      phi = 25, expand.fac = 0.6,
                                      col.colours = NULL, n.col = 100,
                                      main = NULL, ...) {
  .chk_pos(theta1, beta1, theta2, beta2)
  g <- .grid2(n.grid, range)
  G <- as.matrix(expand.grid(g, g))
  z <- relifgmchen(G[, 2], G[, 1], theta1, beta1, theta2, beta2, gamma, omega)
  ttl <- if (is.null(main))
    paste0("IFGM-Chen reliability: gamma=", gamma, ", omega=", omega,
           " (tau=", round(tau_ifgm(gamma, omega), 4), ")") else main
  .surface(z, n.grid, range, theta, phi, expand.fac,
           .palette(col.colours, n.col), ttl, ...)
}

#' @rdname plotbvch
#' @export
plot_contour_bvch <- function(alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5,
                              beta = 1.8, n.grid = 81L, range = 3,
                              n.levels = 12, lab = NULL, add = FALSE,
                              draw.singular = TRUE, main = NULL, ...) {
  .chk_pos(alpha1, alpha2, alpha3, beta)
  g <- .grid2(n.grid, range)
  G <- as.matrix(expand.grid(g, g))
  z <- matrix(dbvch_ac(G[, 2], G[, 1], alpha1, alpha2, alpha3, beta),
              n.grid, n.grid)
  ttl <- if (is.null(main))
    paste0("BvCh ac density: alpha=(", alpha1, ", ", alpha2, ", ", alpha3,
           "), beta=", beta) else main
  if (add) {
    graphics::contour(g, g, z, nlevels = n.levels, add = TRUE,
                      drawlabels = FALSE)
  } else {
    graphics::contour(g, g, z, nlevels = n.levels, xlab = "z1", ylab = "z2",
                      main = if (is.null(lab)) ttl else lab, ...)
  }
  if (draw.singular) {
    graphics::lines(g, g, col = "blue", lwd = 2)
    graphics::legend("topleft", legend = "diagonal (singular support)",
                     col = "blue", lwd = 2, bty = "n")
  }
  invisible(z)
}
#' @noRd
.palette <- function(col.colours, n.col) {
  ## NULL uses the default heat palette; a function is evaluated to n.col
  ## colours; a character vector is interpolated
  if (is.null(col.colours)) return(grDevices::heat.colors(n.col))
  if (is.function(col.colours)) return(col.colours(n.col))
  grDevices::colorRampPalette(col.colours)(n.col)
}

#' @keywords internal
#' @noRd
.surface <- function(z, n.grid, range, theta, phi, expand.fac, pal, main, ...) {
  z <- matrix(z, n.grid, n.grid)
  ## For beta < 1 the density is unbounded at the origin (z^(beta-1)), and for
  ## large beta it underflows to zero on a linear grid.  Clamp so that the
  ## z limits passed to persp are always finite and strictly increasing.
  fin <- z[is.finite(z)]
  if (length(fin) && any(!is.finite(z))) {
    z[!is.finite(z)] <- max(fin) * 1.2
  } else if (!length(fin)) {
    z <- matrix(0, n.grid, n.grid)
  }
  zr <- range(z)
  if (zr[1] >= zr[2]) zr <- c(zr[1] - 0.5, zr[2] + 0.5)
  graphics::persp(x = .grid2(n.grid, range), y = .grid2(n.grid, range), z = z,
                  theta = theta, phi = phi, expand = expand.fac,
                  col = pal, border = NA, xlab = "z1", ylab = "z2",
                  zlab = "", main = main, zlim = zr, ...)
  invisible(z)
}
