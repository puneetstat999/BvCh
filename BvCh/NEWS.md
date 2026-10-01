# BvCh 0.1.0

First release.

## Models

* **BvCh** — Marshall-Olkin type Bivariate Chen distribution. Probability
  density (`dbvch`), survival function (`sbvch`), hazard function (`hbvch`)
  and random generation (`rbvch`), plus the two components of the mixture
  representation (`dbvch_ac`, `dbvch_sing`).
* **BEC** — Bivariate Extended Chen distribution: `dbec`, `sbec`, `hbec`,
  `rbec`, and the componentwise hazard gradients `hbec1` and `hbec2`.
* **IFGM-Chen** — Bivariate Inverse FGM Chen distribution built from the
  Inverse FGM copula and Chen marginals: `relifgmchen`, `difgmchen`,
  `hifgmchen`, `rifgmchen`, together with the copula itself (`cifgm`,
  `cifgmdens`, `rifgm`) and its dependence measures (`tau_ifgm`, `rho_ifgm`).

## Univariate building blocks

* Univariate Chen density, distribution, quantile, survival, hazard,
  hazard-rate-average, conditional reliability and random generation are
  taken from **reliaR** rather than re-implemented.
* The univariate Extended Chen distribution, which reliaR does not provide,
  is implemented here (`sech`, `hech`, `dech`, `pech`, `qech`, `rech`).

## Estimation, comparison and graphics

* Maximum likelihood estimation: `fitbvch`, `fitbec`, `fitifgmchen`, with
  standard errors from the numerical Hessian, and `logLik_bvch`,
  `logLik_bec`, `logLik_ifgmchen` exported for direct use.
* `compare_models` tabulates −log L, AIC and BIC; `ksfit_bvch` performs
  Kolmogorov-Smirnov goodness-of-fit tests for the three marginals.
* Base-graphics surface and contour plots: `plot_surface_bvch`,
  `plot_surface_bec`, `plot_reliability_ifgmchen`, `plot_contour_bvch`.

## Notes

* `BvCh::rchen` draws `reliaR::qchen(runif(n), beta, lambda)` rather than
  calling `reliaR::rchen`, because in reliaR 0.2 that function disagrees
  with the package's own distribution function (Kolmogorov-Smirnov distance
  about 0.38 at `beta = 1.5`, `lambda = 2`, growing to about 0.85 at
  `beta = 2.5`, and correct only when `beta = 1`). See `?chen`.
* `hbvch` implements the hazard as `f/S`; on the diagonal this reproduces
  the source formula `h0(z) = alpha3 * beta * z^(beta-1) * exp(z^beta)`
  exactly.
* With the Inverse FGM copula formula as written in the source document,
  `tau_ifgm` and `rho_ifgm` are non-positive for `gamma, omega >= 0`
  (`tau = -0.0699`, `rho = -0.1049` at `gamma = 0.6`, `omega = 0.4`). The
  source document quotes a positive `tau = 0.156` for the same parameters,
  which is inconsistent with the copula it also quotes. The formula has been
  implemented as given; see `?cifgm`.