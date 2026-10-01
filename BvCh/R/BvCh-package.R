#' Bivariate Chen Distribution: Density, Survival, Hazard, Plots and Estimates
#'
#' Implements probability density, survival and hazard functions, random
#' generation, graphical displays and maximum likelihood estimation for a
#' family of bivariate Chen type lifetime models.
#'
#' @section Models:
#' \describe{
#'   \item{Bivariate Chen (BvCh)}{Marshall-Olkin type model of Gupta, Kumar and
#'     co-authors. With \eqn{U_i \sim Chen(\alpha_i, \beta)} independent and
#'     \eqn{Z_1 = \min(U_1, U_3)}, \eqn{Z_2 = \min(U_2, U_3)}, the pair
#'     \eqn{(Z_1,Z_2)} is \eqn{BvCh(\alpha_1,\alpha_2,\alpha_3,\beta)}. The
#'     distribution is a mixture of an absolutely continuous part of total
#'     mass \eqn{(\alpha_1+\alpha_2)/(\alpha_1+\alpha_2+\alpha_3)} and a
#'     singular part of mass \eqn{\alpha_3/(\alpha_1+\alpha_2+\alpha_3)}
#'     carried on the diagonal.}
#'   \item{Bivariate Extended Chen (BEC)}{Marshall-Olkin type model using
#'     extended Chen latent variables, adding a shape parameter
#'     \eqn{\lambda > 0}.}
#'   \item{Bivariate Inverse FGM Chen (IFGM-Chen)}{Copula construction pairing
#'     the Inverse FGM copula with Chen marginals. This model is absolutely
#'     continuous and has no singular component.}
#' }
#'
#' @section Univariate Chen:
#' All univariate Chen building blocks are taken from the \pkg{reliaR}
#' package rather than being re-implemented. In particular
#' \code{reliaR::schen(x, beta, lambda)} satisfies
#' \code{schen(x, beta, lambda) == exp(lambda * (1 - exp(x^beta))}, which is
#' the \eqn{S_{Ch}(z;\alpha,\beta)} used throughout the bivariate formulae.
#'
#' @docType package
#' @name BvCh-package
#' @aliases BvCh
#' @keywords internal
"_PACKAGE"

## Calls into other packages are made with the explicit `::` operator, but the
## S3 generics used as method targets must be imported so that the
## S3method() directives in NAMESPACE can be resolved.
#' @importFrom stats nobs coef vcov logLik plogis qlogis
NULL
