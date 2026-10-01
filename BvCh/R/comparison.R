#' Model comparison table
#'
#' Collects fitted models into a single table of \eqn{-\log L}, AIC and BIC.
#'
#' @param ... fitted model objects from \code{\link{fitbvch}},
#'   \code{\link{fitbec}} or \code{\link{fitifgmchen}}.
#'
#' @return A data frame with columns \code{Model}, \code{nobs},
#'   \code{neglogL}, \code{AIC} and \code{BIC}, sorted by increasing AIC
#'   so that the best model is on the first row.
#'
#' @examples
#' set.seed(1)
#' a <- fitbvch(rbvch(200, 1.5, 2.5, 0.5, 1.8))
#' b <- fitbec(rbec(200, 1.5, 2.5, 0.5, 1.8, 2))
#' compare_models(a, b)
#'
#' @export
compare_models <- function(...) {
  fits <- list(...)
  if (length(fits) == 0L)
    stop("Supply at least one fitted model.", call. = FALSE)
  nm <- vapply(fits, function(f) {
    if (inherits(f, "fitbvch")) "BvCh"
    else if (inherits(f, "fitbec")) "BEC"
    else if (inherits(f, "fitifgmchen")) "IFGM-Chen"
    else class(f)[1]
  }, character(1))
  df <- data.frame(
    Model = nm,
    nobs = vapply(fits, function(f) f$nobs, numeric(1)),
    neglogL = vapply(fits, function(f) f$negloglik, numeric(1)),
    AIC = vapply(fits, function(f) f$aic, numeric(1)),
    BIC = vapply(fits, function(f) f$bic, numeric(1)),
    stringsAsFactors = FALSE)
  df[order(df$AIC), ]
}

#' Kolmogorov-Smirnov goodness of fit for the BvCh marginals
#'
#' Tests the fitted marginals of a Bivariate Chen fit against the data using
#' \code{\link[reliaR]{ks.chen}}.
#'
#' @details The two components are compared with
#'   \eqn{Chen(\alpha_1+\alpha_3,\beta)} and
#'   \eqn{Chen(\alpha_2+\alpha_3,\beta)}, since \eqn{Z_1=\min(U_1,U_3)} and
#'   \eqn{Z_2=\min(U_2,U_3)}, and their minimum with
#'   \eqn{Chen(\alpha_1+\alpha_2+\alpha_3,\beta)}, which is the marginal of
#'   \eqn{\min(Z_1,Z_2)=\min(U_1,U_2,U_3)}. Each of these is a Chen
#'   distribution with the same shape \eqn{\beta}, because the minimum of
#'   independent Chen variables is Chen in the summed scale parameter.
#'
#' @param fit a fitted model from \code{\link{fitbvch}}.
#' @param x1,x2 numeric vectors of observations, or a two column matrix.
#'
#' @return A data frame with one row per series (\code{Z1}, \code{Z2} and
#'   \code{min(Z1,Z2)}), giving the fitted \code{alpha} and \code{beta}, the
#'   Kolmogorov-Smirnov distance \code{D} and its p-value \code{p.value},
#'   together with the sample mean and standard deviation.
#'
#' @examples
#' set.seed(1)
#' d <- rbvch(200, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
#' fit <- fitbvch(d)
#' ksfit_bvch(fit, d)
#'
#' @export
ksfit_bvch <- function(fit, x1, x2 = NULL) {
  d <- .as_pair(x1, x2, "ksfit_bvch")
  cf <- coef(fit)
  a1 <- cf[["alpha1"]]; a2 <- cf[["alpha2"]]; a3 <- cf[["alpha3"]]
  bt <- cf[["beta"]]
  K <- a1 + a2 + a3
  series <- list(Z1 = d[, 1], Z2 = d[, 2], `min(Z1,Z2)` = pmin(d[, 1], d[, 2]))
  al <- c(a1 + a3, a2 + a3, K)
  out <- do.call(rbind, lapply(seq_along(series), function(i) {
    v <- series[[i]]
    ks <- stats::ks.test(v, reliaR::pchen, bt, al[i])
    data.frame(Series = names(series)[i], alpha = al[i], beta = bt,
               D = unname(ks$statistic), p.value = unname(ks$p.value),
               mean = mean(v), sd = stats::sd(v),
               row.names = NULL, check.names = FALSE)
  }))
  rownames(out) <- NULL
  out
}
