test_that("IFGM copula has the right margins and symmetry", {
  g <- 0.6; w <- 0.4
  u <- seq(0, 1, length.out = 51); v <- seq(0, 1, length.out = 51)
  expect_equal(cifgm(u, 1, g, w), u, tolerance = 1e-14)
  expect_equal(cifgm(1, v, g, w), v, tolerance = 1e-14)
  expect_equal(cifgm(0, v, g, w), rep(0, length(v)), tolerance = 1e-14)
  G <- outer(u, v, function(a, b) cifgm(a, b, g, w))
  expect_equal(G, t(G), tolerance = 1e-14)
  ## independence copula
  expect_equal(cifgm(u, v, 0, 0), u * v, tolerance = 1e-14)
  expect_equal(outer(u, v, function(a, b) cifgm(a, b, 0, 0)),
               outer(u, v, function(a, b) a * b), tolerance = 1e-14)
})

test_that("IFGM copula density equals the second derivative of C", {
  d2 <- function(u, v, g, w, h = 1e-4) {
    (cifgm(u + h, v + h, g, w) - cifgm(u + h, v - h, g, w)
     - cifgm(u - h, v + h, g, w) + cifgm(u - h, v - h, g, w)) / (4 * h * h)
  }
  pts <- expand.grid(u = c(0.2, 0.45, 0.7), v = c(0.15, 0.5, 0.85))
  for (gw in list(c(0, 0), c(0.6, 0.4), c(0.9, 0.1), c(1, 1))) {
    an <- with(pts, cifgmdens(u, v, gw[1], gw[2]))
    nu <- with(pts, d2(u, v, gw[1], gw[2]))
    expect_equal(as.numeric(an), as.numeric(nu), tolerance = 1e-6)
  }
})

test_that("IFGM copula density integrates to one and is positive", {
  for (gw in list(c(0, 0), c(0.6, 0.4), c(0.9, 0.1))) {
    f <- function(u) vapply(u, function(uu) integrate(function(vv)
      cifgmdens(uu, vv, gw[1], gw[2]), 0, 1, rel.tol = 1e-10)$value, numeric(1))
    expect_equal(integrate(f, 0, 1, rel.tol = 1e-10)$value, 1, tolerance = 1e-8)
  }
  g <- seq(0, 1, length.out = 11); v <- seq(0, 1, length.out = 11)
  dens <- outer(g, v, function(a, b) cifgmdens(a, b, 0.6, 0.4))
  expect_true(all(dens > 0))
})

test_that("tau and rho are zero under independence and match simulation", {
  expect_equal(tau_ifgm(0, 0), 0, tolerance = 1e-8)
  expect_equal(rho_ifgm(0, 0), 0, tolerance = 1e-8)
  set.seed(21)
  for (gw in list(c(0.6, 0.4), c(0.2, 0.7))) {
    z <- rifgm(50000, gw[1], gw[2])
    ## Monte Carlo standard errors are about 4*sd(C)/sqrt(n) for tau and
    ## 12*sd(UV)/sqrt(n) for rho, with sd(UV) ~ 1/3
    expect_equal(tau_ifgm(gw[1], gw[2]),
                 4 * mean(cifgm(z[, 1], z[, 2], gw[1], gw[2])) - 1,
                 tolerance = 0.02)
    expect_equal(rho_ifgm(gw[1], gw[2]), 12 * mean(z[, 1] * z[, 2]) - 3,
                 tolerance = 0.06)
  }
})

test_that("rifgm reproduces the copula", {
  set.seed(22)
  z <- rifgm(50000, 0.6, 0.4)
  expect_equal(dim(z), c(50000L, 2L))
  expect_true(all(z >= 0 & z <= 1))
  for (tt in list(c(0.25, 0.25), c(0.5, 0.5), c(0.75, 0.75))) {
    e <- mean(z[, 1] <= tt[1] & z[, 2] <= tt[2])
    expect_equal(e, cifgm(tt[1], tt[2], 0.6, 0.4), tolerance = 0.01)
  }
})

test_that("copula parameters are validated", {
  expect_error(cifgm(0.3, 0.7, -0.1, 0.4))
  expect_error(cifgmdens(0.3, 0.7, 0.6, 1.2))
})

t1 <- 1.2; b1 <- 1.3; t2 <- 1.5; b2 <- 1.1; g <- 0.6; om <- 0.4

test_that("IFGM-Chen density is the mixed derivative of the reliability", {
  d2R <- function(x, y, h = 1e-5) {
    (relifgmchen(x + h, y + h, t1, b1, t2, b2, g, om)
     - relifgmchen(x + h, y - h, t1, b1, t2, b2, g, om)
     - relifgmchen(x - h, y + h, t1, b1, t2, b2, g, om)
     + relifgmchen(x - h, y - h, t1, b1, t2, b2, g, om)) / (4 * h * h)
  }
  for (p in list(c(0.4, 0.7), c(0.8, 0.3), c(1.1, 0.9))) {
    expect_equal(difgmchen(p[1], p[2], t1, b1, t2, b2, g, om), d2R(p[1], p[2]),
                 tolerance = 1e-5)
  }
})

test_that("IFGM-Chen hazard equals f/R", {
  xs <- c(0.4, 0.9, 1.3); ys <- c(0.7, 0.3, 1.1)
  expect_equal(hifgmchen(xs, ys, t1, b1, t2, b2, g, om),
               difgmchen(xs, ys, t1, b1, t2, b2, g, om) /
                 relifgmchen(xs, ys, t1, b1, t2, b2, g, om),
               tolerance = 1e-10)
})

test_that("IFGM-Chen density integrates to one", {
  inner <- function(x) vapply(x, function(zx) integrate(function(zy)
    difgmchen(zx, zy, t1, b1, t2, b2, g, om), 0, 8,
    rel.tol = 1e-10, subdivisions = 3000L)$value, numeric(1))
  tot <- integrate(inner, 0, 8, rel.tol = 1e-8, subdivisions = 3000L)$value
  expect_equal(tot, 1, tolerance = 1e-5)
})

test_that("rifgmchen reproduces the marginals and reliability", {
  set.seed(33)
  s <- rifgmchen(50000, t1, b1, t2, b2, g, om)
  expect_equal(dim(s), c(50000L, 2L))
  q <- c(0.4, 0.7, 1.0)
  expect_lt(max(abs(reliaR::pchen(q, b1, t1) -
    vapply(q, function(v) mean(s[, 1] <= v), numeric(1)))), 0.02)
  expect_lt(max(abs(reliaR::pchen(q, b2, t2) -
    vapply(q, function(v) mean(s[, 2] <= v), numeric(1)))), 0.02)
  pr <- list(c(0.5, 0.8), c(0.9, 0.4))
  emp <- vapply(pr, function(p) mean(s[, 1] > p[1] & s[, 2] > p[2]), numeric(1))
  thS <- vapply(pr, function(p) relifgmchen(p[1], p[2], t1, b1, t2, b2, g, om),
                numeric(1))
  expect_lt(max(abs(thS - emp)), 0.02)
})