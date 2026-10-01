# BvCh

<!-- badges: start -->
<!-- badges: end -->

Probability density, survival function, hazard function, random
generation, graphical displays and maximum likelihood estimation for
**bivariate Chen type lifetime models**.

## Installation

```r
# install.packages("remotes")
remotes::install_github("puneetstat999/BvCh")
```

The univariate Chen building blocks come from
[reliaR](https://cran.r-project.org/package=reliaR), which is installed
automatically.

## The three models

| Model | Construction | Singular part |
|:--|:--|:--|
| **BvCh** | Marshall-Olkin type: $Z_1=\min(U_1,U_3)$, $Z_2=\min(U_2,U_3)$ with $U_i\sim\mathrm{Chen}(\alpha_i,\beta)$ | yes, on the diagonal |
| **BEC** | same construction with extended Chen variables, adding $\lambda$ | yes, on the diagonal |
| **IFGM-Chen** | Inverse FGM copula with Chen marginals | no |

## Quick start

```r
library(BvCh)

## density, survival and hazard
dbvch(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
sbvch(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)
hbvch(0.5, 0.9, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)

## simulation
set.seed(1)
d <- rbvch(500, alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5, beta = 1.8)

## maximum likelihood estimation
fit <- fitbvch(d)
fit
ksfit_bvch(fit, d)

## model comparison
compare_models(fit,
               fitbec(rbec(500, 1.5, 2.5, 0.5, 1.8, 2)),
               fitifgmchen(rifgmchen(500, 1.2, 1.3, 1.5, 1.1, 0.6, 0.4)))

## surfaces
plot_surface_bvch(alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5,
                  beta.vec = c(0.8, 1.8, 5, 10))
plot_reliability_ifgmchen(gamma = 0.6, omega = 0.4)
```

## Univariate Chen

All univariate Chen routines are the [reliaR](https://cran.r-project.org/package=reliaR)
implementations, re-exported for convenience:

```r
dchen(x, beta, lambda)   # density
pchen(q, beta, lambda)   # distribution function
qchen(p, beta, lambda)   # quantile function
rchen(n, beta, lambda)   # random generation
schen(x, beta, lambda)   # survival function
hchen(x, beta, lambda)   # hazard function
```

`reliaR::schen(x, beta, lambda)` equals `exp(lambda * (1 - exp(x^beta)))`,
which is the $S_{Ch}(z;\alpha,\beta)$ used throughout the bivariate
formulae, and `pchen + schen == 1` exactly.

The **extended** Chen distribution, which reliaR does not provide and which
the BEC model needs, is implemented here as `sech`, `hech`, `dech`,
`pech`, `qech`, `rech`.

## Two things worth knowing

1. **`rchen` uses `reliaR::qchen`.** In reliaR 0.2 `reliaR::rchen` does not
   agree with the package's own distribution function — its Kolmogorov-Smirnov
   distance from `pchen` is about 0.38 at `beta = 1.5, lambda = 2` and about
   0.85 at `beta = 2.5`, and it is correct only when `beta = 1`. Since
   `qchen` is the verified analytic inverse of `pchen`, `BvCh::rchen` draws
   `qchen(runif(n), beta, lambda)` so that all six univariate Chen routines
   are mutually consistent. See `?chen`.

2. **The IFGM copula gives negative dependence for positive parameters.**
   With the formula as written in the source document,
   `C(u,v) = uv / sqrt(1 + gamma*(1-u)*(1-v) + omega*u*v*(1-u)*(1-v))`, the
   copula lies below the independence copula `uv`, so `tau_ifgm` and
   `rho_ifgm` are non-positive for `gamma, omega >= 0` (`tau = -0.0699`,
   `rho = -0.1049` at `gamma = 0.6, omega = 0.4`). The source quotes
   `tau = 0.156 > 0` for the same parameters, which contradicts the copula
   formula it also quotes. The formula has been implemented as given; see
   `?cifgm`.

## License

GPL-3

## Authors

* Mukul Bijalwan — <https://orcid.org/0009-0001-3040-6912>
* Puneet Kumar Gupta (maintainer) — <puneetstat999@gmail.com>