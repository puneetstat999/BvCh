a1 <- 1.5; a2 <- 2.5; a3 <- 0.5; bt <- 1.8; lam <- 2
K <- a1 + a2 + a3
x1 <- c(0.4, 0.9, 0.7, 1.3); x2 <- c(0.6, 0.4, 0.7, 1.1)

test_that("BEC survival function matches the piecewise EC product", {
  expect_equal(sbec(x1[1], x2[1], a1, a2, a3, bt, lam),
               sech(x1[1], a1, bt, lam) * sech(x2[1], a2 + a3, bt, lam),
               tolerance = 1e-14)
  expect_equal(sbec(x1[2], x2[2], a1, a2, a3, bt, lam),
               sech(x1[2], a1 + a3, bt, lam) * sech(x2[2], a2, bt, lam),
               tolerance = 1e-14)
  zz <- 0.7
  expect_equal(sbec(zz, zz, a1, a2, a3, bt, lam),
               sech(zz, K, bt, lam), tolerance = 1e-14)
})

test_that("BEC hazard equals f/S", {
  h <- hbec(x1, x2, a1, a2, a3, bt, lam)
  s <- sbec(x1, x2, a1, a2, a3, bt, lam)
  d <- dbec(x1, x2, a1, a2, a3, bt, lam)
  expect_equal(as.numeric(h), as.numeric(d / s), tolerance = 1e-10)
  zz <- 0.7
  expect_equal(hbec(zz, zz, a1, a2, a3, bt, lam),
               a3 * bt * lam * zz^(bt - 1) * exp(zz^bt) /
                 (1 + lam * expm1(zz^bt)), tolerance = 1e-10)
})

test_that("BEC hazard gradients multiply to the joint hazard off diagonal", {
  off <- x1 != x2
  h12 <- hbec1(x1, x2, a1, a2, a3, bt, lam) * hbec2(x1, x2, a1, a2, a3, bt, lam)
  expect_equal(h12[off], hbec(x1, x2, a1, a2, a3, bt, lam)[off],
               tolerance = 1e-10)
})

test_that("BEC region masses equal alpha_i / (sum alpha) and sum to one", {
  m1 <- integrate(function(z) vapply(z, function(v)
    integrate(function(v2) dbec(v, v2, a1, a2, a3, bt, lam), v, Inf,
              rel.tol = 1e-10, subdivisions = 2000L)$value, numeric(1)),
    0, Inf, rel.tol = 1e-8, subdivisions = 2000L)$value
  m2 <- integrate(function(z) vapply(z, function(v) if (v > 0)
    integrate(function(v2) dbec(v, v2, a1, a2, a3, bt, lam), 0, v,
              rel.tol = 1e-10, subdivisions = 2000L)$value else 0, numeric(1)),
    0, Inf, rel.tol = 1e-8, subdivisions = 2000L)$value
  u <- seq(1e-8, 8, length.out = 200001L)
  m3 <- sum(dbec(u, u, a1, a2, a3, bt, lam)) * (u[2] - u[1])
  expect_equal(m1, a1 / K, tolerance = 1e-5)
  expect_equal(m2, a2 / K, tolerance = 1e-5)
  expect_equal(m3, a3 / K, tolerance = 1e-5)
  expect_equal(m1 + m2 + m3, 1, tolerance = 1e-5)
})

test_that("rbec reproduces the marginals and joint survival", {
  set.seed(11)
  sim <- rbec(100000, a1, a2, a3, bt, lam)
  expect_equal(dim(sim), c(100000L, 2L))
  q <- c(0.4, 0.6, 0.8, 1.0)
  th <- 1 - sech(q, a1, bt, lam) * sech(q, a3, bt, lam)
  em <- vapply(q, function(v) mean(sim[, 1] <= v), numeric(1))
  expect_lt(max(abs(th - em)), 0.01)
  pr <- list(c(0.5, 0.8), c(0.9, 0.4))
  emp <- vapply(pr, function(p) mean(sim[, 1] > p[1] & sim[, 2] > p[2]), numeric(1))
  thS <- vapply(pr, function(p) sech(p[1], a1, bt, lam) * sech(p[2], a2, bt, lam) *
    sech(max(p[1], p[2]), a3, bt, lam), numeric(1))
  expect_lt(max(abs(thS - emp)), 0.01)
})