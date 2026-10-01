a1 <- 1.5; a2 <- 2.5; a3 <- 0.5; bt <- 1.8
K <- a1 + a2 + a3
z1 <- c(0.4, 0.9, 0.7, 1.3); z2 <- c(0.6, 0.4, 0.7, 1.1)

test_that("BvCh survival function matches the product form", {
  S <- sbvch(z1, z2, a1, a2, a3, bt)
  Sref <- reliaR::schen(z1, bt, a1) * reliaR::schen(z2, bt, a2) *
    reliaR::schen(pmax(z1, z2), bt, a3)
  expect_equal(as.numeric(S), as.numeric(Sref), tolerance = 1e-14)
  expect_true(all(S >= 0 & S <= 1))
})

test_that("BvCh density branches agree with the paper's f1, f2, f0", {
  ## z1 < z2 -> f1
  f1 <- a1 * (a2 + a3) * bt^2 * z1[1]^(bt - 1) * z2[1]^(bt - 1) *
    exp(z1[1]^bt + z2[1]^bt + a1 * (1 - exp(z1[1]^bt)) +
      (a2 + a3) * (1 - exp(z2[1]^bt)))
  expect_equal(dbvch(z1[1], z2[1], a1, a2, a3, bt), f1, tolerance = 1e-12)
  ## z1 > z2 -> f2
  f2 <- a2 * (a1 + a3) * bt^2 * z1[2]^(bt - 1) * z2[2]^(bt - 1) *
    exp(z1[2]^bt + z2[2]^bt + a2 * (1 - exp(z2[2]^bt)) +
      (a1 + a3) * (1 - exp(z1[2]^bt)))
  expect_equal(dbvch(z1[2], z2[2], a1, a2, a3, bt), f2, tolerance = 1e-12)
  ## z1 = z2 -> f0
  zz <- 0.77
  f0 <- a3 * bt * zz^(bt - 1) *
    exp(zz^bt + K * (1 - exp(zz^bt)))
  expect_equal(dbvch(zz, zz, a1, a2, a3, bt), f0, tolerance = 1e-12)
})

test_that("BvCh hazard equals f/S and matches h0 on the diagonal", {
  h <- hbvch(z1, z2, a1, a2, a3, bt)
  s <- sbvch(z1, z2, a1, a2, a3, bt)
  d <- dbvch(z1, z2, a1, a2, a3, bt)
  expect_equal(as.numeric(h), as.numeric(d / s), tolerance = 1e-10)
  zz <- 0.77
  expect_equal(hbvch(zz, zz, a1, a2, a3, bt), a3 * bt * zz^(bt - 1) * exp(zz^bt),
               tolerance = 1e-12)
})

test_that("Proposition 1 mixture representation holds", {
  ac <- dbvch_ac(z1, z2, a1, a2, a3, bt)
  sg <- dbvch_sing(z1, z2, a1, a2, a3, bt)
  mix <- (a1 + a2) / K * ac + a3 / K * sg
  expect_equal(as.numeric(mix),
               as.numeric(dbvch(z1, z2, a1, a2, a3, bt)), tolerance = 1e-14)
  ## the ac part vanishes on the diagonal, the singular part only lives there
  expect_equal(unname(ac[z1 == z2]), 0)
  expect_gt(sg[z1 == z2][1], 0)
  expect_equal(unname(sg[z1 != z2]), c(0, 0, 0))
})

test_that("BvCh density integrates to one over the three regions", {
  ## mass of the region z1 < z2, of z1 > z2 and of the diagonal
  m1 <- integrate(function(z) vapply(z, function(v)
    integrate(function(v2) dbvch(v, v2, a1, a2, a3, bt), v, Inf,
              rel.tol = 1e-10, subdivisions = 2000L)$value, numeric(1)),
    0, Inf, rel.tol = 1e-8, subdivisions = 2000L)$value
  m2 <- integrate(function(z) vapply(z, function(v) if (v > 0)
    integrate(function(v2) dbvch(v, v2, a1, a2, a3, bt), 0, v,
              rel.tol = 1e-10, subdivisions = 2000L)$value else 0, numeric(1)),
    0, Inf, rel.tol = 1e-8, subdivisions = 2000L)$value
  u <- seq(1e-8, 6, length.out = 200001L)
  m0 <- a3 / K * sum(reliaR::dchen(u, bt, K)) * (u[2] - u[1])
  expect_equal(m1, a1 / K, tolerance = 1e-5)
  expect_equal(m2, a2 / K, tolerance = 1e-5)
  expect_equal(m0, a3 / K, tolerance = 1e-4)
  expect_equal(m1 + m2 + m0, 1, tolerance = 1e-4)
})

test_that("rbvch reproduces the theoretical marginals and tie frequency", {
  set.seed(42)
  sim <- rbvch(100000, a1, a2, a3, bt)
  expect_equal(dim(sim), c(100000L, 2L))
  q <- seq(0.4, 1.6, length.out = 6)
  th <- 1 - reliaR::schen(q, bt, a1) * reliaR::schen(q, bt, a3)
  em <- vapply(q, function(v) mean(sim[, 1] <= v), numeric(1))
  expect_lt(max(abs(th - em)), 0.01)
  th2 <- 1 - reliaR::schen(q, bt, a2) * reliaR::schen(q, bt, a3)
  em2 <- vapply(q, function(v) mean(sim[, 2] <= v), numeric(1))
  expect_lt(max(abs(th2 - em2)), 0.01)
  expect_equal(mean(sim[, 1] == sim[, 2]), a3 / K, tolerance = 0.01)
})

test_that("vectorisation and recycling work", {
  v <- dbvch(c(0.4, 0.9), 0.6, a1, a2, a3, bt)
  expect_length(v, 2)
  expect_equal(v[1], dbvch(0.4, 0.6, a1, a2, a3, bt))
  expect_equal(v[2], dbvch(0.9, 0.6, a1, a2, a3, bt))
  expect_equal(dbvch(-1, 2, a1, a2, a3, bt), 0)
})