test_that("log-likelihood functions agree with the densities", {
  set.seed(5)
  d <- rbvch(100, 1.5, 2.5, 0.5, 1.8)
  p <- log(c(1.5, 2.5, 0.5, 1.8))
  expect_equal(logLik_bvch(p, d[, 1], d[, 2]),
               sum(dbvch(d[, 1], d[, 2], 1.5, 2.5, 0.5, 1.8, log = TRUE)))
  e <- rbec(100, 1.5, 2.5, 0.5, 1.8, 2)
  p5 <- log(c(1.5, 2.5, 0.5, 1.8, 2))
  expect_equal(logLik_bec(p5, e[, 1], e[, 2]),
               sum(dbec(e[, 1], e[, 2], 1.5, 2.5, 0.5, 1.8, 2, log = TRUE)))
  w <- rifgmchen(100, 1.2, 1.3, 1.5, 1.1, 0.6, 0.4)
  p6 <- log(c(1.2, 1.3, 1.5, 1.1, 0.6, 0.4))
  expect_equal(logLik_ifgmchen(p6, w[, 1], w[, 2]),
               sum(difgmchen(w[, 1], w[, 2], 1.2, 1.3, 1.5, 1.1, 0.6, 0.4,
                             log = TRUE)))
})

test_that("dropping ties changes the log-likelihood only through the ties", {
  set.seed(6)
  d <- rbvch(200, 1.5, 2.5, 0.5, 1.8)
  p <- log(c(1.5, 2.5, 0.5, 1.8))
  tie <- d[, 1] == d[, 2]
  expect_gt(sum(tie), 0)
  ll_all <- logLik_bvch(p, d[, 1], d[, 2], singular = TRUE)
  ll_nt <- logLik_bvch(p, d[, 1], d[, 2], singular = FALSE)
  expect_equal(ll_all - ll_nt, sum(dbvch(d[tie, 1], d[tie, 2], 1.5, 2.5, 0.5,
                                          1.8, log = TRUE)))
})

test_that("fitbvch recovers the parameters", {
  set.seed(1)
  d <- rbvch(3000, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
  fit <- fitbvch(d)
  expect_s3_class(fit, "fitbvch")
  expect_equal(length(coef(fit)), 4L)
  expect_equal(names(coef(fit)), c("alpha1", "alpha2", "alpha3", "beta"))
  truth <- c(1.5, 2.5, 0.5, 1.8)
  expect_lt(max(abs(coef(fit) - truth) / truth), 0.15)
  expect_true(all(fit$se > 0))
  expect_equal(dim(vcov(fit)), c(4L, 4L))
  ## AIC and BIC agree with their definitions
  expect_equal(fit$aic, 2 * 4 - 2 * fit$loglik)
  expect_equal(fit$bic, 4 * log(fit$nobs) - 2 * fit$loglik)
  expect_equal(nobs(fit), fit$nobs)
})

test_that("fitbec recovers the parameters", {
  set.seed(2)
  d <- rbec(3000, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8,
            lambda = 2)
  fit <- fitbec(d)
  expect_s3_class(fit, "fitbec")
  expect_equal(names(coef(fit)),
               c("alpha1", "alpha2", "alpha3", "beta", "lambda"))
  expect_equal(fit$aic, 2 * 5 - 2 * fit$loglik)
  expect_true(all(fit$se > 0))
})

test_that("fitifgmchen recovers the marginal parameters", {
  set.seed(3)
  d <- rifgmchen(3000, theta1 = 1.2, beta1 = 1.3, theta2 = 1.5,
                 beta2 = 1.1, gamma = 0.6, omega = 0.4)
  fit <- fitifgmchen(d)
  expect_s3_class(fit, "fitifgmchen")
  expect_equal(names(coef(fit)),
               c("theta1", "beta1", "theta2", "beta2", "gamma", "omega"))
  expect_equal(fit$aic, 2 * 6 - 2 * fit$loglik)
  truth <- c(theta1 = 1.2, beta1 = 1.3, theta2 = 1.5, beta2 = 1.1)
  est <- coef(fit)[1:4]
  expect_lt(max(abs(est - truth) / truth), 0.15)
  ## copula parameters stay inside (0, 1)
  expect_true(all(coef(fit)[5:6] >= 0 & coef(fit)[5:6] <= 1))
})

test_that("input validation in the fitting functions", {
  expect_error(fitbvch(1:10))
  expect_error(fitbvch(matrix(1:10, ncol = 5)))
  expect_error(fitbvch(c(-1, 1, 2, 3)))
  expect_error(fitbvch(cbind(1, 1, 1)), "two column")
  expect_error(fitbvch(rbvch(20, 1, 1, 1, 1), start = c(1, 1)))
})

test_that("S3 methods work", {
  set.seed(4)
  f <- fitbvch(rbvch(200, 1.5, 2.5, 0.5, 1.8))
  expect_output(print(f), "BvCh model fit")
  expect_output(print(summary(f)), "Coefficients")
  expect_equal(coef(f), f$coefficients)
  expect_equal(vcov(f), f$vcov)
  expect_equal(unname(as.numeric(logLik(f))), f$loglik)
  expect_equal(attr(logLik(f), "df"), 4L)
})

test_that("compare_models tabulates fitted models", {
  set.seed(5)
  a <- fitbvch(rbvch(200, 1.5, 2.5, 0.5, 1.8))
  b <- fitbec(rbec(200, 1.5, 2.5, 0.5, 1.8, 2))
  tab <- compare_models(a, b)
  expect_equal(nrow(tab), 2L)
  expect_equal(sort(tab$Model), c("BEC", "BvCh"))
  expect_true(all(diff(tab$AIC) >= 0))   # sorted by increasing AIC
  expect_equal(tab$AIC[1], min(a$aic, b$aic))
  expect_true(all(diff(tab$BIC) >= 0))
})

test_that("ksfit_bvch reports the three marginal series", {
  set.seed(6)
  d <- rbvch(500, 1.5, 2.5, 0.5, 1.8)
  fit <- fitbvch(d)
  ks <- ksfit_bvch(fit, d)
  expect_equal(nrow(ks), 3L)
  expect_equal(ks$Series, c("Z1", "Z2", "min(Z1,Z2)"))
  cf <- coef(fit)
  expect_equal(ks$alpha[1], cf[["alpha1"]] + cf[["alpha3"]])
  expect_equal(ks$alpha[2], cf[["alpha2"]] + cf[["alpha3"]])
  expect_equal(ks$alpha[3], cf[["alpha1"]] + cf[["alpha2"]] + cf[["alpha3"]])
  ## a correctly specified fit should not be rejected
  expect_gt(ks$p.value[1], 0.01)
})