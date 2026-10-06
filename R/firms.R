#' Simulated firm panel with an omitted variable
#'
#' A simulated data set of 250 firms observed for 20 years, used in the
#' examples and the vignette. It is a variant of the simulated example in
#' Lonati and Wulff (2026, Table 1). The outcome \code{y} depends on the
#' predictor of interest \code{x}, whose true effect is 0.25, and on two
#' firm-level characteristics of equal strength: \code{c}, which is observed,
#' and \code{q}, which is meant to be omitted from the analysis. Because
#' \code{q} also drives \code{x}, a regression of \code{y} on \code{x} and
#' \code{c} overstates the effect of \code{x}. The within-firm components of
#' \code{x} and of the disturbance are persistent (first-order autoregressive
#' with coefficient 0.9), so standard errors have to be clustered by firm.
#'
#' The variable \code{q} is included so that every sensitivity result can be
#' checked against the regression that contains it. "An omitted variable as
#' strong as \code{c}" (\code{benchmark_covariates = "c"}, \code{kd = 1}) is,
#' by construction, a variable like \code{q}.
#'
#' @format A data frame with 5,000 rows and 6 columns:
#' \describe{
#'   \item{\code{firm}}{firm identifier, 1 to 250 (the cluster)}
#'   \item{\code{year}}{year within firm, 1 to 20}
#'   \item{\code{y}}{outcome}
#'   \item{\code{x}}{predictor of interest; its true effect on \code{y} is 0.25}
#'   \item{\code{c}}{observed firm-level control}
#'   \item{\code{q}}{unobserved firm-level variable, omitted in the analysis}
#' }
#' @details
#' The data were generated as follows (the script is \code{data-raw/firms.R}
#' in the package repository; the Stata command ships the same data as
#' \code{firms.dta}): \code{q} and \code{c} are independent standard normal
#' draws per firm; \code{v} and \code{u} are independent first-order
#' autoregressive series within each firm with persistence 0.9 and stationary
#' variance 1; \code{x = 0.5 q + 0.5 c + v} and
#' \code{y = 0.25 x + 0.5 q + 0.5 c + u}. The draws of \code{q}, \code{v} and
#' \code{u} are those of seed 116 of the paper's script. Compared with the
#' paper's Table 1, the control \code{c} is added so that the benchmark
#' approach applies, and the coefficients of \code{c} and \code{q} are 0.5
#' instead of 1 and 0.7.
#'
#' @references
#' Lonati, S., & Wulff, J. N. (2026). Why you should not use the ITCV with
#' robust standard errors (and what to do instead). \emph{Academy of
#' Management Proceedings}, 2026(1).
#' \doi{10.5465/AMPROC.2026.247bp}
#'
#' @examples
#' data("firms", package = "bootmakr")
#' head(firms)
#'
#' # The regression that omits q overstates the effect of x (true value 0.25);
#' # the one that includes q recovers it.
#' coef(lm(y ~ x + c, data = firms))["x"]
#' coef(lm(y ~ x + c + q, data = firms))["x"]
"firms"
