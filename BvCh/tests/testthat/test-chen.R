test_that("univariate Chen is consistent with reliaR", {
  x <- c(0.2, 0.5, 1.0, 2.0)
  b <- 1.5; la <- 2
  expect_equal(dchen(x, b, la), reliaR::dchen(x, b, la))
  expect_equal(schen(x, b, la), reliaR::schen(x, b, la))
  expect_equal(pchen(x, b, la), reliaR::pchen(x, b, la))
  ## schen is exactly exp(lambda * (1 - exp(x^beta))), the S_Ch of the papers
  expect_equal(schen(x, b, la), exp(la * (1 - exp(x^b))), tolerance = 1e-14)
  expect_equal(pchen(x, b, la) + schen(x, b, la), rep(1, length(x)),
               tolerance = 1e-14)
  expect_equal(hchen(x, b, la), dchen(x, b, la) / schen(x, b, la),
               tolerance = 1e-14)
  expect_equal(pchen(qchen(c(0.1, 0.5, 0.9), b, la), b, la),
               c(0.1, 0.5, 0.9), tolerance = 1e-12)
})

test_that("rchen follows the Chen distribution", {
  ## reliaR::rchen itself is inconsistent with its own CDF (see the docs);
  ## BvCh::rchen must be the corrected version.
  set.seed(1)
  s <- rchen(5000, 1.5, 2)
  ps <- sort(pchen(s, 1.5, 2))
  n <- length(ps); i <- seq_len(n)
  ks <- max(max(i / n - ps), max(ps - (i - 1) / n))
  expect_lt(ks, 0.03)
  expect_gt(stats::ks.test(s, reliaR::pchen, 1.5, 2)$p.value, 0.01)
})

test_that("univariate Extended Chen is internally consistent", {
  x <- c(0.2, 0.5, 1.0)
  p <- c(0.05, 0.3, 0.6, 0.95)
  expect_equal(pech(x, 1.5, 1.5, 2) + sech(x, 1.5, 1.5, 2), rep(1, 3))
  expect_equal(dech(x, 1.5, 1.5, 2), sech(x, 1.5, 1.5, 2) * hech(x, 1.5, 1.5, 2))
  expect_equal(pech(qech(p, 1.5, 1.5, 2), 1.5, 1.5, 2), p, tolerance = 1e-12)
  expect_equal(qchen <- pchen(qech(p, 1.5, 1.5, 2), 1.5, 1.5, 2), pchen(qech(p, 1.5, 1.5, 2), 1.5, 1.5, 2))
  set.seed(1)
  s <- rech(5000, 1.5, 1.5, 2)
  ps <- sort(pech(s, 1.5, 1.5, 2))
  n <- length(ps); i <- seq_len(n)
  expect_lt(max(max(i / n - ps), max(ps - (i - 1) / n)), 0.03)
})

test_that("parameters must be positive", {
  expect_error(dbvch(1, 2, -1, 1, 1, 1))
  expect_error(sbvch(1, 2, 1, 1, 1, 0))
  expect_error(dbec(1, 2, 1, 1, 1, 1, -1))
  expect_error(sech(1, 1, 1, 0))
})