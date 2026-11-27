#' Bootstrap Inference for Sensemakr Sensitivity Analysis
#'
#' Wraps \code{\link[sensemakr]{sensemakr}} (Cinelli & Hazlett, 2020) in a
#' bootstrap loop. For an omitted variable \code{kd} times as strong as an
#' observed benchmark covariate, \code{bootmakr()} computes the bias-adjusted
#' estimate of the treatment coefficient and obtains its standard error,
#' confidence interval and p-value by simple, stratified or cluster bootstrap,
#' so that inference does not rest on the conventional OLS variance formula.
#'
#' @details
#' In every replication the data are resampled (rows, whole clusters, or
#' within strata), the regression is refitted, and the benchmark bounds and
#' the bias-adjusted estimate are recomputed. The bootstrap distribution
#' therefore reflects sampling uncertainty in the regression coefficient
#' \emph{and} in the estimated strength of the benchmark. The standard error
#' is the standard deviation of that distribution, the confidence interval
#' its percentile interval, and the p-value is twice the smaller of the
#' shares of replications at or below zero and at or above zero.
#'
#' With a length-two \code{cluster}, three bootstraps are combined with the
#' Cameron, Gelbach & Miller (2011) subtractive variance estimator, and the
#' confidence interval and p-value use a normal approximation.
#'
#' Below the table of adjusted estimates, the print method reports the
#' \strong{benchmark strength}: the partial R2 of the benchmark with the
#' treatment (given the other covariates) and with the outcome (given the
#' treatment and the covariates), and sensemakr's implied bound on the
#' omitted variable at each \code{kd}. For the omitted variable it also
#' reports its \strong{impact}, the product of its partial correlations with
#' the outcome and with the treatment, both given the covariates only. This
#' is the scale on which the impact threshold of a confounding variable
#' (ITCV) is stated, so a reader used to the ITCV can judge how strong the
#' assumed omitted variable is. The partial correlation with the outcome
#' given the covariates is recovered from sensemakr's bound (which conditions
#' on the treatment as well) with the recursion formula for partial
#' correlations, with the signs set so that the omitted variable biases the
#' estimate away from zero (towards zero when \code{reduce = FALSE}). All of
#' these are functions of the data alone. No standard error enters them, so
#' they do not depend on the variance estimator, and none of them is a
#' threshold. They are stored in \code{$benchmark_strength}.
#'
#' The print method then reports the \strong{breakdown kd}: the smallest
#' \code{kd}, with \code{ky = kd}, at which the estimate is no longer
#' significant at \code{alpha}, for the benchmark the table of adjusted
#' estimates refers to. Two numbers are given. The first is judged by the
#' confidence interval that \code{sensemakr::ovb_bounds()} reports, which
#' rests on the conventional OLS standard error. The second is judged by the
#' bootstrap test of the table: running \code{bootmakr()} with the same
#' \code{seed} and a \code{kd} marginally below that value returns a p-value
#' below \code{alpha}, and marginally above it a p-value of at least
#' \code{alpha}. It is computed without refitting, from the coefficient, its
#' scale and the benchmark's two partial R2 in every replicate
#' (\code{$boot_base} and \code{$boot_bench}), and it does not depend on the
#' \code{kd} of the call. sensemakr's bound exists only up to a largest
#' \code{kd}, at which the implied partial R2 of the omitted variable with the
#' treatment reaches one. That limit is printed. It differs between
#' replicates, and a replicate beyond its own limit is left out, as it is in
#' the table; a note says how many are left out at the reported value. A
#' benchmark that is strongly related to the treatment has a low limit, and
#' close to the limit the bias grows without bound, so a breakdown kd near the
#' limit reflects the limit more than the estimate. With \code{weights} both
#' numbers use the weighted partial R2 of the benchmark, and with
#' \code{semiweights} the bound of Wainstein and Hazlett (2026), whose limit
#' is \code{(1 - R2w) / R2s}. The breakdown kd is
#' available for \code{reduce = TRUE} and when the table reports a benchmark
#' rather than a manual R2. It is stored in \code{$breakdown}.
#'
#' The print method also reports the \strong{robustness value} of the
#' treatment coefficient, which needs no benchmark: the strength an omitted
#' variable must have, as a partial R2 with the treatment and with the
#' outcome alike, to overturn the result. Three numbers are given. The
#' first two are \code{sensemakr}'s: the strength that brings the adjusted
#' estimate to zero, and the strength at which the adjusted estimate is no
#' longer significant at \code{alpha} when judged with the conventional OLS
#' standard error. The third is the bootstrap counterpart of the second:
#' the smallest strength at which the bootstrap test used for the table of
#' adjusted estimates stops rejecting.
#' Running \code{bootmakr()} with the same \code{seed} and with \code{r2dz.x}
#' and \code{r2yz.dx} both set marginally above that value returns a p-value
#' of at least \code{alpha}; marginally below it, a p-value below
#' \code{alpha}. As in the table, the omitted variable is taken to pull the
#' coefficient of each replicate towards zero, whatever its sign in that
#' replicate. When many replicates have a coefficient of the opposite sign to
#' the estimate, as with a weak effect, this convention raises the bootstrap
#' value, and a bootstrap value above \code{sensemakr}'s should then not be
#' read as greater robustness.
#' It is
#' computed from the coefficient and the residual scale of every replicate,
#' without refitting, uses every replicate in which the regression could be
#' fitted, and is available for \code{q = 1} and \code{reduce = TRUE}. It
#' carries the same notion of uncertainty as the table: the sampling
#' variability of the adjusted estimate for an omitted variable of given
#' strength. It is not a statement about the robust standard error that the
#' regression including the omitted variable would have had. The three
#' numbers are stored in \code{$robustness}.
#'
#' The breakdown kd and the bootstrap robustness value are set by the tail of
#' the bootstrap distribution: with \code{alpha = 0.05}, by the 2.5 percent of
#' replicates closest to zero. They therefore need more replications than a
#' standard error does. Use 1,000 or more, and expect them still to differ by
#' a few percent from one seed to the next. The sign convention described for
#' the robustness value applies to the breakdown kd as well, so with a weak
#' effect the bootstrap value can exceed the conventional one for that reason
#' alone.
#'
#' After the robustness values the print method sets the \strong{standard
#' errors of the adjusted estimate} side by side, one row per \code{kd}. The
#' first is the bootstrap standard error of the table. The second, shown when
#' a benchmark is re-estimated in every replicate, is the bootstrap standard
#' error with the strength of the omitted variable held at its value in the
#' original fit. It is computed from the same replicates, and the difference
#' between the two is the part of the uncertainty that comes from estimating
#' the benchmark. The third is \code{sensemakr}'s adjusted standard error,
#' which rests on the conventional OLS standard error. The fourth is the
#' closed form of Wainstein and Hazlett (2026): the robust standard error of
#' the original fit multiplied by \code{sqrt((1 - r2yz.dx) / (1 - r2dz.x))},
#' with one more regressor in its finite-sample correction. The robust
#' standard error is HC1 without clusters, CR1 with one cluster variable and,
#' as an extension made here, the Cameron, Gelbach and Miller combination of
#' CR1 variances with two; strata are ignored. The finite-sample corrections
#' of the variance are \code{n / (n - k)} for HC1 and
#' \code{G / (G - 1) * (n - 1) / (n - k)} for CR1, with \code{n} the number
#' of observations with a positive weight, \code{G} the number of clusters
#' and \code{k} the number of regressors, which for the adjusted standard
#' error includes the omitted variable (equation 58 of Wainstein and Hazlett,
#' 2026). The closed form assumes that
#' the omitted variable shrinks every residual of the outcome and of the
#' treatment by the same proportion, and Acedanski (2026) shows that an
#' omitted variable can move a robust standard error by much more than that.
#' The last two standard errors refer to the regression that would include the
#' omitted variable. The bootstrap refers to the sampling variability of the
#' adjusted estimate, so the columns answer different questions and need not
#' agree. They are stored in \code{$se_compare}.
#'
#' The print method ends with the \strong{ratio of the robust to the
#' conventional standard error} of the treatment coefficient in the unadjusted
#' regression, and with the change in that ratio when each benchmark covariate
#' is taken out of the regression and put back (Acedanski, 2026, Section 6.2).
#' A ratio away from one means that robust and conventional inference already
#' differ before any omitted variable is considered. The change is a yardstick
#' for how far one more covariate might move the ratio. Both are stored in
#' \code{$se_ratio}, whose \code{cap} is the ratio plus the largest positive
#' change.
#'
#' The bootstrap procedure follows Cinelli, Ferwerda and Hazlett (2024,
#' Appendix C); Lonati and Wulff (2026) show why it is needed when the
#' regression calls for robust standard errors.
#'
#' @param formula A formula for the OLS regression: outcome on the left,
#'   treatment and covariates on the right.
#' @param data A data frame.
#' @param treat Character: name of the treatment variable (the coefficient
#'   whose sensitivity is assessed).
#' @param benchmark_covariates Character vector: observed covariate(s) used as
#'   benchmark for the strength of the omitted variable. If several are
#'   supplied, the adjusted estimates refer to the first one.
#' @param gbenchmark_covariates Character vector: a group of covariates used
#'   jointly as benchmark. Uses \code{group_partial_r2} to compute joint R2
#'   values.
#' @param kd,ky Numeric vectors of benchmark multipliers: the omitted variable
#'   is assumed \code{kd} times as strong as the benchmark in explaining the
#'   treatment and \code{ky} times as strong in explaining the outcome.
#'   \code{kd} defaults to 1 (exactly as strong as the benchmark) and
#'   \code{ky} to \code{kd}. A vector of several values gives one row of
#'   results per value (a "kd sweep").
#' @param q,alpha Numeric: proportion of the effect to be explained away
#'   (passed to sensemakr) and significance level of the bootstrap confidence
#'   interval (default 0.05).
#' @param r2dz.x,r2yz.dx Optional manual partial R2 values of the omitted
#'   variable with the treatment and with the outcome (instead of a
#'   benchmark).
#' @param bound_label Character label for bounds table.
#' @param reduce Logical: should the omitted variable reduce the absolute
#'   value of the estimate (default \code{TRUE})?
#' @param bounds_row Integer: which row of sensemakr's bounds table to use
#'   when a single \code{kd} is supplied (default 1).
#' @param reps Integer: bootstrap replications (default 1000).
#' @param seed Integer or \code{NULL}. No seed is set unless one is supplied,
#'   so results vary from run to run without it.
#' @param cluster Column name(s) or vector for the cluster bootstrap.
#'   Length 1 (e.g. \code{"firm"}) gives a one-way cluster bootstrap with
#'   percentile CIs. Length 2 (e.g. \code{c("firm", "year")}) triggers a
#'   two-way cluster bootstrap using the Cameron, Gelbach & Miller (2011)
#'   subtractive variance estimator: three bootstraps are run (on each
#'   dimension and on their intersection) and combined as
#'   \code{V = V_G + V_H - V_(G n H)}. CIs and p-values for the two-way case
#'   use a normal approximation \code{estimate +/- z * SE}, or a t distribution
#'   with \code{twoway_crit = "t"}. Two-way clustering
#'   is incompatible with \code{strata}.
#' @param strata,weights Column name or vector: strata for a stratified
#'   bootstrap, and regression weights. Weights are treated as fixed: each
#'   unit is resampled together with its weight, and the weights are not
#'   re-estimated within replications. They enter the outcome model and the
#'   treatment model that measures the benchmark's strength, so every
#'   partial R2 is a weighted one.
#'   \code{weights} may instead be a function that takes a data frame and
#'   returns one non-negative weight per row, for example inverse propensity
#'   score weights. It is called on the data and again on every resample, so
#'   that the weights are re-estimated within replications (Option 1 in
#'   Wainstein & Hazlett, 2026). A resample on which the function fails counts
#'   as a failed replication. Unless \code{semiweights} is given, a benchmark
#'   is measured with the weights themselves; see \code{semiweights} for a
#'   benchmark covariate that enters the weights. A weight function together
#'   with a benchmark and no \code{semiweights} gives a warning.
#' @param dots Deprecated; use \code{progress} instead.
#' @param progress Logical: show a progress bar (default TRUE).
#' @param converge \code{TRUE}, \code{FALSE}, or
#'   \code{list(minreps, stepsize, threshold)}: an informal convergence check
#'   that recomputes the standard error and p-value on the first
#'   \code{minreps}, \code{minreps + stepsize}, ... replications and
#'   summarises how much they still move, overall and from \code{threshold}
#'   replications onwards.
#' @param verbose Logical; currently unused.
#' @param fix_benchmark Logical. If \code{TRUE}, the partial R2 of the
#'   benchmark with the treatment and with the outcome are taken from the
#'   original fit in every replicate instead of being re-estimated there, so
#'   that the omitted variable has the same strength throughout. The table,
#'   the breakdown kd and \code{$boot_bench} then refer to that fixed
#'   strength. This is the bootstrap of Wainstein and Hazlett (2026), who hold
#'   the two partial R2 of the omitted variable fixed. It leaves out the
#'   uncertainty about how strong the benchmark is. Default \code{FALSE}.
#' @param twoway_crit With two-way clustering, the reference distribution of
#'   the confidence interval, the p-value, the breakdown kd and the bootstrap
#'   robustness value: \code{"normal"} (the default, as in earlier versions)
#'   or \code{"t"}, a t distribution with \code{min(G, H) - 1} degrees of
#'   freedom, where \code{G} and \code{H} are the numbers of clusters in the
#'   two dimensions (Cameron, Gelbach & Miller, 2011). Ignored otherwise.
#' @param semiweights Optional semi-weights for the benchmark (Wainstein &
#'   Hazlett, 2026): weights formed in the same way as \code{weights} but
#'   without the benchmark covariate. Weights that balance a covariate
#'   between treated and control units remove its association with the
#'   treatment, so that its weighted partial R2 with the treatment is close
#'   to zero and says little about how strongly a variable like it can
#'   confound. With \code{semiweights} that partial R2 is measured twice,
#'   with the weights (\code{R2w}) and with the semi-weights (\code{R2s}),
#'   and the bound on the omitted variable is that of Wainstein and Hazlett
#'   (their equations 25 and 26), which starts from
#'   \code{r2dz.x = kd * R2s / (1 - R2w)}. \code{kd} then compares the
#'   omitted variable under the weights with the benchmark under the
#'   semi-weights. This is their definition and not that of Cinelli and
#'   Hazlett (2020); the two coincide, and the bound is sensemakr's, when the
#'   semi-weights equal the weights. The outcome side is measured with the
#'   weights throughout. Supply the semi-weights in the same form as
#'   \code{weights}: a column name or a vector, resampled with the units,
#'   when \code{weights} is one; a function of the data, called again on
#'   every resample, when \code{weights} is a function. With several
#'   \code{benchmark_covariates} give a list with one entry per covariate,
#'   each formed without that covariate; \code{gbenchmark_covariates} takes
#'   one entry, formed without the whole group. A resample on which a
#'   semi-weight function fails counts as a failed replication. Default
#'   \code{NULL}: the benchmark is measured with the weights. With a weight
#'   function this default gives a warning, because estimated weights usually
#'   depend on the covariates; if yours do not depend on the benchmark, pass
#'   the weight function as \code{semiweights} too, which gives the same
#'   results without the warning.
#' @return An object of class \code{"bootmakr"}, a list with (among others)
#'   \describe{
#'     \item{\code{results}}{data frame with one row per \code{kd}: adjusted
#'       \code{estimate}, bootstrap \code{se}, \code{ci_lower},
#'       \code{ci_upper}, \code{pvalue}, the implied strength of the
#'       omitted variable (\code{r2dz.x}, \code{r2yz.dx}) the row is based
#'       on, and its \code{impact}.}
#'     \item{\code{benchmark_strength}}{list with \code{$observed} (partial
#'       R2 of the benchmark with treatment and outcome, and their square
#'       roots; with \code{semiweights} also \code{r2dxj.x_semi}, the
#'       partial R2 with the treatment under the semi-weights),
#'       \code{$implied} (sensemakr's bound on the omitted variable
#'       at each \code{kd}/\code{ky}: \code{r2dz.x}, \code{r2yz.dx}, their
#'       square roots \code{r_dz.x}, \code{r_yz.dx}, the partial correlation
#'       with the outcome given the covariates only, \code{r_yz.x}, and the
#'       \code{impact}, \code{r_yz.x * r_dz.x}) and \code{$r_yd.x} (the
#'       partial correlation of the outcome with the treatment given the
#'       covariates, used in that recursion).}
#'     \item{\code{boot_samples}}{matrix of bootstrap draws of the adjusted
#'       estimate, one column per \code{kd} (a list of three such matrices
#'       with two-way clustering).}
#'     \item{\code{breakdown}}{list with the breakdown kd of the benchmark
#'       the table refers to: \code{kd} (bootstrap), \code{kd_conventional}
#'       (conventional standard error), \code{kd_limit} (the largest kd the
#'       benchmark allows in the original fit), \code{benchmark},
#'       \code{n_used} and \code{n_reps} (replicates with a bound at the
#'       reported value, and replicates with a benchmark at all), and a
#'       \code{note} where one applies.}
#'     \item{\code{boot_bench}}{matrix with the partial R2 of that benchmark
#'       with the treatment (\code{r2dxj.x}) and with the outcome
#'       (\code{r2yxj.dx}) in every replicate (a list of three such matrices
#'       with two-way clustering); \code{NA} without a benchmark. With
#'       \code{semiweights} a third column, \code{r2dxj.x_semi}, holds the
#'       partial R2 with the treatment under the semi-weights.}
#'     \item{\code{robustness}}{list with the robustness values of the
#'       treatment coefficient: \code{rv_q} and \code{rv_qa} from
#'       \code{sensemakr} (the latter at \code{alpha}, with the
#'       conventional standard error) and \code{rv_qa_boot}, the bootstrap
#'       robustness value; also \code{n_used}, the number of replicates it
#'       is based on, and a \code{note} when it is not available.}
#'     \item{\code{boot_base}}{matrix with the unadjusted treatment
#'       coefficient (\code{estimate}) and its \code{scale}, the standard
#'       error times the square root of the residual degrees of freedom, in
#'       every replicate (a list of three such matrices with two-way
#'       clustering). The adjusted estimate for any strength of the omitted
#'       variable can be computed from these two columns.}
#'     \item{\code{se_compare}}{list with \code{$table}, one row per \code{kd}
#'       with the standard errors of the adjusted estimate (\code{bootstrap},
#'       \code{bootstrap_fixed}, \code{conventional}, \code{robust}),
#'       \code{$unadjusted} (conventional and robust standard error of the
#'       treatment coefficient), \code{$robust_type} and a \code{note} where
#'       one applies.}
#'     \item{\code{se_ratio}}{list with \code{omega}, the ratio of the robust
#'       to the conventional standard error of the treatment coefficient,
#'       \code{$benchmarks} (the ratio without each benchmark and the
#'       \code{change} when it is restored) and \code{cap}.}
#'     \item{\code{convergence}}{convergence diagnostics, if requested.}
#'     \item{\code{N}, \code{N_reps}, \code{N_successful}, \code{N_fail},
#'       \code{N_clust}}{sample size, replications requested, completed and
#'       failed, and number of clusters.}
#'     \item{\code{sensemakr_orig}}{the \code{sensemakr} object of the
#'       original fit. With \code{weights} its benchmark bounds can differ
#'       from those in \code{results}: \code{sensemakr} (checked in version
#'       0.1.5) measures the benchmark's association with the treatment
#'       without the weights.}
#'   }
#'   Use \code{print()} for the results table and \code{plot()} for the
#'   kd-sweep, convergence and histogram plots.
#'
#' @references
#' Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2011). Robust inference
#' with multiway clustering. \emph{Journal of Business & Economic
#' Statistics}, 29(2), 238-249.
#'
#' Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity: Extending
#' omitted variable bias. \emph{Journal of the Royal Statistical Society:
#' Series B}, 82(1), 39-67.
#'
#' Cinelli, C., Ferwerda, J., & Hazlett, C. (2024). sensemakr: Sensitivity
#' analysis tools for OLS in R and Stata. \emph{Observational Studies},
#' 10(2), 93-127.
#'
#' Lonati, S., & Wulff, J. N. (2026). Why you should not use the ITCV with
#' robust standard errors (and what to do instead). \emph{Academy of
#' Management Proceedings}, 2026(1).
#' \doi{10.5465/AMPROC.2026.247bp}
#'
#' Wainstein, L., & Hazlett, C. (2026). Sensitivity of weighted least squares
#' estimators to omitted variables. arXiv:2508.02954.
#'
#' Acedanski, J. (2026). Robustness values with sandwich standard errors.
#' Working paper, SSRN 7639842.
#'
#' @seealso \code{\link{print.bootmakr}}, \code{\link{plot.bootmakr}},
#'   \code{\link{firms}}, \code{\link[sensemakr]{sensemakr}}
#'
#' @examples
#' # Simulated firm panel: the true effect of x on y is 0.25, and a firm-level
#' # variable q, omitted here, biases the regression of y on x and c (see
#' # ?firms). An omitted variable as strong as the observed control c (kd = 1,
#' # the default); firms are resampled as whole clusters. 200 replications
#' # keep the example quick; use 1,000 or more in practice.
#' data("firms", package = "bootmakr")
#' out <- bootmakr(y ~ x + c, data = firms, treat = "x",
#'                 benchmark_covariates = "c", cluster = "firm",
#'                 reps = 200, seed = 1, progress = FALSE)
#' out
#'
#' # The strength of the benchmark and of the omitted variable it implies,
#' # with the impact of the latter on the scale of the ITCV
#' out$benchmark_strength
#'
#' # Several strengths at once: where does the effect stop being significant?
#' sweep <- bootmakr(y ~ x + c, data = firms, treat = "x",
#'                   benchmark_covariates = "c", kd = c(0.5, 1, 1.5, 2),
#'                   cluster = "firm", reps = 200, seed = 1, progress = FALSE)
#' sweep$results
#' plot(sweep, type = "kd_sweep")
#'
#' @importFrom sensemakr sensemakr group_partial_r2 ovb_partial_r2_bound adjusted_estimate partial_r2 ovb_bounds
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline hist legend mtext par plot points segments
#' @importFrom stats coef complete.cases lm model.matrix nobs pnorm qnorm quantile reformulate sd terms var
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
bootmakr <- function(formula,
                     data,
                     treat,
                     benchmark_covariates  = NULL,
                     gbenchmark_covariates = NULL,
                     kd = 1,
                     ky = NULL,
                     q = 1,
                     alpha = 0.05,
                     r2dz.x = NULL,
                     r2yz.dx = NULL,
                     bound_label = NULL,
                     reduce = TRUE,
                     bounds_row = 1,
                     reps = 1000,
                     seed = NULL,
                     cluster = NULL,
                     strata = NULL,
                     weights = NULL,
                     progress = TRUE,
                     dots = 0,
                     converge = FALSE,
                     verbose = FALSE,
                     fix_benchmark = FALSE,
                     twoway_crit = c("normal", "t"),
                     semiweights = NULL) {

  cl <- match.call()
  if (is.null(ky)) ky <- kd
  n_kd <- length(kd)
  twoway_crit <- match.arg(twoway_crit)
  if (!is.logical(fix_benchmark) || length(fix_benchmark) != 1 || is.na(fix_benchmark))
    stop("`fix_benchmark` must be TRUE or FALSE.")

  # Validate: need exactly one of benchmark / gbenchmark / manual R2
  n_bench_args <- sum(!is.null(benchmark_covariates),
                      !is.null(gbenchmark_covariates),
                      !is.null(r2dz.x) && !is.null(r2yz.dx))
  if (n_bench_args == 0)
    stop("Supply one of: benchmark_covariates, gbenchmark_covariates, or r2dz.x + r2yz.dx")
  if (n_bench_args > 1)
    warning("Multiple benchmark sources supplied; gbenchmark_covariates takes priority.")

  use_gbench <- !is.null(gbenchmark_covariates)
  bench_label <- if (use_gbench) {
    paste(gbenchmark_covariates, collapse = ", ")
  } else if (!is.null(benchmark_covariates)) {
    paste(benchmark_covariates, collapse = ", ")
  } else {
    bound_label %||% "manual"
  }

  # Resolve cluster / strata / weights
  cluster_list  <- .resolve_cluster(cluster, data)
  strata_vec    <- .resolve_var(strata, data, "strata")
  # `weights` may be a function of the data: it is then called on the data
  # and again on every resample, so the weights are re-estimated within
  # replications instead of being resampled with their units
  weight_fun    <- if (is.function(weights)) weights else NULL
  weight_vec    <- if (is.null(weight_fun)) .resolve_var(weights, data, "weights")
                   else .call_weight_fun(weight_fun, data)
  # Semi-weights of the benchmark(s): the weights formed without the
  # benchmark. One entry per benchmark covariate (one for a grouped
  # benchmark), of the same kind as `weights`: vectors are resampled with
  # their units, functions are called again on every resample.
  semi          <- .resolve_semiweights(semiweights, data, weights, weight_fun,
                                        benchmark_covariates, gbenchmark_covariates)
  semi_orig     <- .semi_values(semi, data)
  # Weights estimated from the data usually come from the covariates, and
  # then they balance the benchmark between treated and control units: its
  # partial R2 with the treatment under the weights is close to zero and the
  # bound built on it understates the bias. Columns and vectors of weights
  # cannot be told apart in this respect (survey weights are fine), so the
  # warning is given for weight functions only.
  if (!is.null(weight_fun) && is.null(semi) &&
      (!is.null(benchmark_covariates) || !is.null(gbenchmark_covariates)))
    warning("`weights` is a function and no `semiweights` are given, so the benchmark is ",
            "measured with the weights themselves. Weights built from the covariates ",
            "balance the benchmark between treated and control units; its partial R2 with ",
            "the treatment is then close to zero and the adjusted estimate understates the ",
            "bias. Supply `semiweights`, the weights formed without the benchmark. If the ",
            "weights do not depend on the benchmark, pass the weight function as ",
            "`semiweights` as well: the results are the same and this warning is not given.",
            call. = FALSE)
  n_cluster_dims <- if (is.null(cluster_list)) 0L else length(cluster_list)
  is_twoway     <- n_cluster_dims == 2L
  if (is_twoway && !is.null(strata_vec))
    stop("`strata` is not supported with two-way clustering.")

  conv_opts <- .parse_converge(converge, reps)
  if (!is.null(seed)) set.seed(seed)

  # ---- Fit original model ----
  fit_orig <- .fit_lm(formula, data, weight_vec)
  N <- nobs(fit_orig)

  # ---- Original point estimates ----
  obs_estimates <- .get_adjusted_estimates(
    fit_orig, data, formula, treat, weight_vec,
    benchmark_covariates, gbenchmark_covariates,
    kd, ky, q, alpha, r2dz.x, r2yz.dx, bound_label, reduce, bounds_row,
    semi = semi_orig
  )

  # Also store the full sensemakr object for the first kd (for reference)
  sm_orig <- .make_sensemakr_orig(
    fit_orig, treat, benchmark_covariates, gbenchmark_covariates,
    kd, ky, q, alpha, r2dz.x, r2yz.dx, bound_label, reduce
  )

  # ---- Benchmark strength (descriptive; no standard errors involved) ----
  # Partial R2 of the benchmark with treatment and outcome, and sensemakr's
  # implied omitted-variable strength at each kd/ky. Computed once on the
  # full sample; not part of the bootstrap loop and draws no random numbers.
  bench_strength <- .benchmark_strength(
    fit_orig, data, formula, treat, weight_vec,
    benchmark_covariates, gbenchmark_covariates,
    kd, ky, r2dz.x, r2yz.dx, bound_label, reduce, semi = semi_orig
  )

  # ---- Set up resampling ----
  if (is_twoway) {
    intersection_vec <- paste(cluster_list[[1]], cluster_list[[2]], sep = "_._")
    resample_info_G  <- .setup_resampling(data, cluster_list[[1]],  NULL)
    resample_info_H  <- .setup_resampling(data, cluster_list[[2]],  NULL)
    resample_info_GH <- .setup_resampling(data, intersection_vec,   NULL)
    .warn_few_clusters(resample_info_G$n_clust,
                       if (is.character(cluster)) cluster[1] else "dim 1")
    .warn_few_clusters(resample_info_H$n_clust,
                       if (is.character(cluster)) cluster[2] else "dim 2")
  } else {
    single_clust  <- if (n_cluster_dims == 1L) cluster_list[[1]] else NULL
    resample_info <- .setup_resampling(data, single_clust, strata_vec)
  }

  # Which benchmark the table of adjusted estimates refers to (NA if it
  # reports a manual R2 row or there is no benchmark)
  bench_idx <- .table_benchmark_index(benchmark_covariates, gbenchmark_covariates,
                                      r2dz.x, n_kd, bounds_row)

  # With fix_benchmark the partial R2s of the benchmark(s) in the original
  # fit are used in every replicate instead of being re-estimated there
  base_fixed <- if (isTRUE(fix_benchmark))
    .bench_base(fit_orig, data, formula, treat, weight_vec,
                benchmark_covariates, gbenchmark_covariates, semi_orig) else NULL

  # Reference distribution of the two-way test: normal, or t with
  # min(G, H) - 1 degrees of freedom
  crit_df <- if (is_twoway && twoway_crit == "t")
    min(resample_info_G$n_clust, resample_info_H$n_clust) - 1 else Inf
  if (is_twoway && is.finite(crit_df) && crit_df < 1)
    stop("`twoway_crit = \"t\"` needs at least two clusters in each dimension.")

  # ---- Bootstrap loop(s) ----
  loop_args <- list(
    reps = reps, data = data, formula = formula, treat = treat,
    weight_vec = weight_vec, weight_fun = weight_fun, base_fixed = base_fixed,
    semi = semi,
    benchmark_covariates = benchmark_covariates,
    gbenchmark_covariates = gbenchmark_covariates,
    kd = kd, ky = ky, q = q, alpha = alpha,
    r2dz.x = r2dz.x, r2yz.dx = r2yz.dx,
    bound_label = bound_label, reduce = reduce, bounds_row = bounds_row,
    bench_idx = bench_idx,
    show_progress = isTRUE(progress)
  )

  if (is_twoway) {
    nm1 <- if (is.character(cluster)) cluster[1] else "dim 1"
    nm2 <- if (is.character(cluster)) cluster[2] else "dim 2"
    out_G  <- do.call(.run_bootstrap_loop,
                      c(list(resample_info = resample_info_G,
                             label = sprintf("[1/3] cluster: %s", nm1)),
                        loop_args))
    out_H  <- do.call(.run_bootstrap_loop,
                      c(list(resample_info = resample_info_H,
                             label = sprintf("[2/3] cluster: %s", nm2)),
                        loop_args))
    out_GH <- do.call(.run_bootstrap_loop,
                      c(list(resample_info = resample_info_GH,
                             label = sprintf("[3/3] intersection (%s x %s)", nm1, nm2)),
                        loop_args))

    boot_samples <- list(G = out_G$boot_mat, H = out_H$boot_mat, GH = out_GH$boot_mat)
    boot_base    <- list(G = out_G$base_mat, H = out_H$base_mat, GH = out_GH$base_mat)
    boot_bench   <- list(G = out_G$bench_mat, H = out_H$bench_mat, GH = out_GH$bench_mat)
    results      <- .compute_twoway_stats(boot_samples, obs_estimates, kd, alpha, crit_df)
    n_fail       <- out_G$n_fail + out_H$n_fail + out_GH$n_fail
    n_successful <- min(sum(complete.cases(boot_samples$G)),
                        sum(complete.cases(boot_samples$H)),
                        sum(complete.cases(boot_samples$GH)))
    n_clust      <- c(G  = resample_info_G$n_clust,
                      H  = resample_info_H$n_clust,
                      GH = resample_info_GH$n_clust)
    method       <- "cgm_twoway"
  } else {
    out <- do.call(.run_bootstrap_loop,
                   c(list(resample_info = resample_info, label = NULL), loop_args))
    boot_samples <- out$boot_mat
    boot_base    <- out$base_mat
    boot_bench   <- out$bench_mat
    results      <- .compute_boot_stats(boot_samples, obs_estimates, kd, alpha)
    n_fail       <- out$n_fail
    n_successful <- sum(complete.cases(boot_samples))
    n_clust      <- resample_info$n_clust
    method       <- "percentile"
  }
  results <- .attach_r2(results, bench_strength, n_kd, r2dz.x, r2yz.dx,
                        bounds_row, reduce)

  # ---- Robustness value (does not depend on the benchmark or on kd) ----
  robustness <- .robustness(fit_orig, treat, q, alpha, reduce, boot_base, is_twoway,
                            crit_df)

  # ---- Breakdown kd: where the bootstrap test stops rejecting ----
  breakdown <- .breakdown(fit_orig, treat, alpha, reduce, kd, ky, bench_idx,
                          bench_strength, boot_base, boot_bench, is_twoway,
                          crit_df)

  # ---- Standard errors side by side, and the robust-to-conventional ratio ----
  # Descriptive companions of the table: none of them draws random numbers
  # or changes anything above.
  se_groups <- .se_groups(fit_orig, cluster_list, nrow(data))
  bench_reestimated <- !isTRUE(fix_benchmark) && !is.na(bench_idx) &&
    bench_strength$type %in% c("single", "group")
  se_compare <- tryCatch(
    .se_compare(fit_orig, treat, results, se_groups, boot_base, reduce,
                is_twoway, bench_reestimated, !is.null(strata_vec),
                if (is.character(cluster)) cluster else NULL),
    error = function(e) NULL)
  se_ratio <- tryCatch(
    .se_ratio(fit_orig, treat, se_groups, benchmark_covariates,
              gbenchmark_covariates, if (is.character(cluster)) cluster else NULL),
    error = function(e) NULL)

  conv_out <- NULL
  if (conv_opts$do_converge) {
    conv_out <- if (is_twoway) {
      .convergence_diagnostics_twoway(boot_samples, obs_estimates[1], alpha, conv_opts,
                                      crit_df)
    } else {
      .convergence_diagnostics(boot_samples[, 1], obs_estimates[1], conv_opts)
    }
  }

  structure(
    list(
      results       = results,
      boot_samples  = boot_samples,
      boot_base     = boot_base,
      boot_bench    = boot_bench,
      robustness    = robustness,
      breakdown     = breakdown,
      se_compare    = se_compare,
      se_ratio      = se_ratio,
      fix_benchmark = isTRUE(fix_benchmark),
      weights_reestimated = !is.null(weight_fun),
      semiweights   = !is.null(semi),
      twoway_crit   = if (is_twoway) twoway_crit else NULL,
      crit_df       = crit_df,
      reduce        = reduce,
      convergence   = conv_out,
      call          = cl,
      method        = method,
      N             = N,
      N_reps        = reps,
      N_successful  = n_successful,
      N_fail        = n_fail,
      N_clust       = n_clust,
      cluster_names = if (is.character(cluster)) cluster else NULL,
      kd            = kd,
      ky            = ky,
      alpha         = alpha,
      treat         = treat,
      benchmark_covariates  = benchmark_covariates,
      gbenchmark_covariates = gbenchmark_covariates,
      bench_label   = bench_label,
      benchmark_strength = bench_strength,
      sensemakr_orig = sm_orig
    ),
    class = "bootmakr"
  )
}


# ==============================================================================
# Core extraction: unified pathway for regular & grouped benchmarks
# ==============================================================================

#' Get adjusted estimates for all kd values from a fitted model
#' @noRd
.get_adjusted_estimates <- function(fit, data, formula, treat, weight_vec,
                                    benchmark_covariates, gbenchmark_covariates,
                                    kd, ky, q, alpha, r2dz.x, r2yz.dx,
                                    bound_label, reduce, bounds_row,
                                    base = NULL, semi = NULL) {
  # `base`: the benchmark's partial R2s if the caller has already computed
  # them for this fit with .bench_base() (the bootstrap loop does, to store
  # them); NULL computes them here.
  # `semi`: semi-weights of the rows of `data`, one vector per benchmark, or
  # NULL; see .ovb_bound().
  n_kd <- length(kd)

  if (!is.null(gbenchmark_covariates)) {
    # ---- Grouped benchmark pathway ----
    # 1.-2. Partial R2 of Y with Z_group given D, X (outcome model) and of
    #       D with Z_group given X (treatment model) -- see .group_base_r2()
    if (is.null(base))
      base <- .group_base_r2(fit, data, formula, treat, weight_vec, gbenchmark_covariates,
                             semi)
    r2yxj_base <- base$r2yxj.dx
    r2dxj_base <- base$r2dxj.x

    # 3. For each kd/ky, use ovb_partial_r2_bound for proper nonlinear scaling,
    #    then compute adjusted estimate with the scaled R2 values
    coefs  <- summary(fit)$coefficients
    est    <- coefs[treat, "Estimate"]
    se     <- coefs[treat, "Std. Error"]
    dof    <- fit$df.residual

    vals <- vapply(seq_len(n_kd), function(i) {
      bounds <- .ovb_bound(
        r2dxj.x = r2dxj_base, r2yxj.dx = r2yxj_base,
        kd = kd[i], ky = ky[i], bound_label = "group",
        r2dxj.x_semi = base$r2dxj.x_semi
      )
      sensemakr::adjusted_estimate(
        estimate = est, se = se, dof = dof,
        r2dz.x  = bounds$r2dz.x,
        r2yz.dx = bounds$r2yz.dx,
        reduce  = reduce
      )
    }, numeric(1))
    return(vals)
  }

  # ---- Regular benchmark / manual R2 pathway ----
  # Computed step by step instead of through sensemakr::sensemakr().
  # sensemakr's benchmarking routine refits the treatment model as an
  # UNWEIGHTED lm(d ~ XX + 0) whatever weights `fit` carries, so with
  # `weights` the treatment-side partial R2 of the benchmark, and the bound
  # kd * R2 / (1 - R2) built on it, ignored the weights while the outcome
  # side used them (and the grouped pathway above used them on both sides).
  # Here the treatment model carries the weights of `fit`; see
  # .treatment_model_sm(). For an unweighted `fit` every step is the call
  # sensemakr itself makes, so unweighted results are unchanged. Row order
  # mirrors sensemakr's bounds table: manual bound(s) first, then the
  # benchmarks, kd varying within benchmark.
  coefs <- summary(fit)$coefficients
  if (!(treat %in% rownames(coefs)))
    stop("Variables not found in model: ", treat)
  est   <- coefs[treat, "Estimate"]
  se    <- coefs[treat, "Std. Error"]
  dof   <- fit$df.residual
  vals  <- numeric(0)

  if (!is.null(r2dz.x)) {
    r2yz_manual <- if (is.null(r2yz.dx)) r2dz.x else r2yz.dx
    if (any(r2dz.x < 0 | r2dz.x > 1 | r2yz_manual < 0 | r2yz_manual > 1))
      stop("partial R2 must be a number between zero and one")
    vals <- c(vals, sensemakr::adjusted_estimate(
      estimate = est, se = se, dof = dof,
      r2dz.x = r2dz.x, r2yz.dx = r2yz_manual, reduce = reduce))
  }

  if (!is.null(benchmark_covariates)) {
    if (is.null(base)) base <- .single_base_r2(fit, treat, benchmark_covariates, semi)
    for (j in seq_along(benchmark_covariates)) {
      bnd <- .ovb_bound(
        r2dxj.x = base$r2dxj.x[j], r2yxj.dx = base$r2yxj.dx[j],
        kd = kd, ky = ky, bound_label = benchmark_covariates[j],
        r2dxj.x_semi = base$r2dxj.x_semi[j])
      vals <- c(vals, sensemakr::adjusted_estimate(
        estimate = est, se = se, dof = dof,
        r2dz.x = bnd$r2dz.x, r2yz.dx = bnd$r2yz.dx, reduce = reduce))
    }
  }

  if (length(vals) == 0) return(rep(NA_real_, n_kd))
  if (n_kd > 1) vals[seq_len(n_kd)] else vals[bounds_row]
}


#' Build a sensemakr object for the original fit (for reference / printing)
#' @noRd
.make_sensemakr_orig <- function(fit, treat, benchmark_covariates,
                                 gbenchmark_covariates, kd, ky, q, alpha,
                                 r2dz.x, r2yz.dx, bound_label, reduce) {
  args <- list(model = fit, treatment = treat)
  if (!is.null(benchmark_covariates)) args$benchmark_covariates <- benchmark_covariates
  args$kd <- kd; args$ky <- ky; args$q <- q; args$alpha <- alpha
  if (!is.null(r2dz.x))     args$r2dz.x     <- r2dz.x
  if (!is.null(r2yz.dx))    args$r2yz.dx     <- r2yz.dx
  if (!is.null(bound_label)) args$bound_label <- bound_label
  args$reduce <- reduce
  tryCatch(do.call(sensemakr::sensemakr, args), error = function(e) NULL)
}


# ==============================================================================
# Benchmark strength (descriptive)
# ==============================================================================

#' Group partial R2s of a set of benchmark covariates: with the outcome given
#' treatment and the other covariates (r2yxj.dx, from the outcome model) and
#' with the treatment given the other covariates (r2dxj.x, from a treatment
#' model treat ~ all other RHS terms). Shared by the bootstrap pathway and the
#' descriptive block so that both report the same numbers. With semi-weights
#' (`semi`, a list holding one vector for the group) the treatment model is
#' fitted a second time with them, and r2dxj.x_semi is the group's partial R2
#' in that fit.
#' @noRd
.group_base_r2 <- function(fit, data, formula, treat, weight_vec, gbenchmark_covariates,
                           semi = NULL) {
  r2yxj.dx <- sensemakr::group_partial_r2(fit, covariates = gbenchmark_covariates)
  rhs_vars     <- attr(terms(formula), "term.labels")
  rhs_no_treat <- setdiff(rhs_vars, treat)
  treat_formula <- reformulate(rhs_no_treat, response = treat)
  fit_d <- .fit_lm(treat_formula, data, weight_vec)
  r2dxj.x <- sensemakr::group_partial_r2(fit_d, covariates = gbenchmark_covariates)
  out <- list(r2dxj.x = unname(r2dxj.x), r2yxj.dx = unname(r2yxj.dx))
  if (!is.null(semi)) {
    fit_s <- .fit_lm(treat_formula, data, semi[[1]])
    out$r2dxj.x_semi <- unname(
      sensemakr::group_partial_r2(fit_s, covariates = gbenchmark_covariates))
  }
  out
}

#' Treatment model as sensemakr's benchmarking routine builds it (model matrix
#' of the outcome model, treatment regressed on all other columns, no
#' intercept beyond the model-matrix constant), except that it carries the
#' weights of the outcome model. sensemakr fits it unweighted even when the
#' outcome model is weighted. `fit$weights` has one entry per model-matrix
#' row, so the two stay aligned when rows are dropped for missing values.
#' `w` replaces those weights, for the treatment model under semi-weights;
#' it must have one entry per model-matrix row (see .rows_used()).
#' @noRd
.treatment_model_sm <- function(fit, treat, w = fit$weights) {
  m  <- model.matrix(fit)
  d  <- m[, treat]
  XX <- m[, !(colnames(m) %in% treat), drop = FALSE]
  w_fit <- w
  if (is.null(w_fit)) lm(d ~ XX + 0) else lm(d ~ XX + 0, weights = w_fit)
}

#' Entries of a vector with one value per row of the data that belong to the
#' rows lm() kept in `fit`.
#' @noRd
.rows_used <- function(v, fit) {
  dropped <- as.integer(fit$na.action)
  if (length(dropped)) v[-dropped] else v
}

#' Partial R2s of each single benchmark covariate: with the outcome given the
#' treatment and the other covariates (r2yxj.dx, from the outcome model) and
#' with the treatment given the other covariates (r2dxj.x, from the treatment
#' model, weighted like the outcome model). Shared by the bootstrap pathway
#' and the descriptive block so that both report the same numbers. With
#' semi-weights (`semi`, a list with one vector per benchmark covariate) the
#' treatment model is fitted once more per benchmark with that benchmark's
#' semi-weights, and r2dxj.x_semi is its partial R2 in that fit.
#' @noRd
.single_base_r2 <- function(fit, treat, benchmark_covariates, semi = NULL) {
  r2yxj.dx <- sensemakr::partial_r2(fit, covariates = benchmark_covariates)
  tm       <- .treatment_model_sm(fit, treat)
  r2dxj.x  <- sensemakr::partial_r2(tm, covariates = paste0("XX", benchmark_covariates))
  out <- list(r2dxj.x = unname(r2dxj.x), r2yxj.dx = unname(r2yxj.dx))
  if (!is.null(semi)) {
    out$r2dxj.x_semi <- vapply(seq_along(benchmark_covariates), function(j) {
      tm_s <- .treatment_model_sm(fit, treat, w = .rows_used(semi[[j]], fit))
      unname(sensemakr::partial_r2(tm_s, covariates = paste0("XX", benchmark_covariates[j])))
    }, numeric(1))
  }
  out
}

#' Descriptive benchmark strength.
#'
#' Returns a list with
#'   $type      "single" (benchmark_covariates), "group" (gbenchmark_covariates),
#'              "manual" (r2dz.x + r2yz.dx supplied) or "unavailable"
#'   $observed  data.frame, one row per benchmark covariate (or one row for the
#'              group): r2dxj.x, r2yxj.dx and their square roots r_dxj.x,
#'              r_yxj.dx (absolute partial correlations). NULL for "manual".
#'   $implied   data.frame, one row per benchmark x kd: kd, ky, r2dz.x, r2yz.dx
#'              (sensemakr's bound on the omitted variable's strength), their
#'              square roots r_dz.x, r_yz.dx, the partial correlation with the
#'              outcome given the covariates only, r_yz.x, and the impact
#'              r_yz.x * r_dz.x (see .impact()).
#'   $r_yd.x    signed partial correlation of the outcome with the treatment
#'              given the covariates (from the fitted model).
#' All quantities are functions of the data alone (no standard errors), so
#' they are unaffected by heteroskedasticity- or cluster-robust inference.
#' Nothing here is a t-implied correlation or a threshold.
#' @noRd
.benchmark_strength <- function(fit, data, formula, treat, weight_vec,
                                benchmark_covariates, gbenchmark_covariates,
                                kd, ky, r2dz.x, r2yz.dx, bound_label,
                                reduce = TRUE, semi = NULL) {
  n_kd <- length(kd)
  out  <- tryCatch({
    if (!is.null(gbenchmark_covariates)) {
      base  <- .group_base_r2(fit, data, formula, treat, weight_vec, gbenchmark_covariates,
                              semi)
      label <- paste(gbenchmark_covariates, collapse = " + ")
      observed <- data.frame(benchmark = label,
                             r2dxj.x = base$r2dxj.x, r2yxj.dx = base$r2yxj.dx,
                             stringsAsFactors = FALSE)
      observed$r2dxj.x_semi <- base$r2dxj.x_semi
      implied <- do.call(rbind, lapply(seq_len(n_kd), function(i) {
        b <- .ovb_bound(
          r2dxj.x = base$r2dxj.x, r2yxj.dx = base$r2yxj.dx,
          kd = kd[i], ky = ky[i], bound_label = "group",
          r2dxj.x_semi = base$r2dxj.x_semi
        )
        data.frame(benchmark = label, kd = kd[i], ky = ky[i],
                   r2dz.x = b$r2dz.x, r2yz.dx = b$r2yz.dx,
                   stringsAsFactors = FALSE)
      }))
      list(type = "group", observed = observed, implied = implied)
    } else if (!is.null(benchmark_covariates)) {
      # Observed strength, computed as in sensemakr's benchmarking routine
      # but with the treatment model weighted like the outcome model
      base <- .single_base_r2(fit, treat, benchmark_covariates, semi)
      observed <- data.frame(benchmark = benchmark_covariates,
                             r2dxj.x = base$r2dxj.x, r2yxj.dx = base$r2yxj.dx,
                             stringsAsFactors = FALSE)
      observed$r2dxj.x_semi <- base$r2dxj.x_semi
      # Implied strength of the omitted variable: sensemakr's bound formula
      # applied to those partial R2s (rows ordered benchmark-by-benchmark, kd
      # within benchmark -- the same order used for the adjusted estimates).
      # Not sensemakr::ovb_bounds(fit, ...), which would refit the treatment
      # model without the weights.
      implied <- do.call(rbind, lapply(seq_along(benchmark_covariates), function(j) {
        b <- .ovb_bound(
          r2dxj.x = base$r2dxj.x[j], r2yxj.dx = base$r2yxj.dx[j],
          kd = kd, ky = ky, bound_label = benchmark_covariates[j],
          r2dxj.x_semi = base$r2dxj.x_semi[j]
        )
        data.frame(benchmark = benchmark_covariates[j], kd = kd, ky = ky,
                   r2dz.x = b$r2dz.x, r2yz.dx = b$r2yz.dx,
                   stringsAsFactors = FALSE)
      }))
      list(type = "single", observed = observed, implied = implied)
    } else {
      implied <- data.frame(benchmark = bound_label %||% "manual",
                            kd = NA_real_, ky = NA_real_,
                            r2dz.x = r2dz.x, r2yz.dx = r2yz.dx,
                            stringsAsFactors = FALSE)
      list(type = "manual", observed = NULL, implied = implied)
    }
  }, error = function(e) list(type = "unavailable", observed = NULL, implied = NULL,
                              error = conditionMessage(e)))

  # Partial correlation of the outcome with the treatment given the
  # covariates (signed), needed to put the omitted variable on the ITCV scale.
  out$r_yd.x <- tryCatch({
    b <- coef(fit)[[treat]]
    unname(sign(b) * sqrt(sensemakr::partial_r2(fit, covariates = treat)))
  }, error = function(e) NA_real_)

  if (!is.null(out$implied)) {
    out$implied$r_dz.x  <- sqrt(out$implied$r2dz.x)
    out$implied$r_yz.dx <- sqrt(out$implied$r2yz.dx)
    imp <- .impact(out$implied$r2dz.x, out$implied$r2yz.dx, out$r_yd.x, reduce)
    out$implied$r_yz.x  <- imp$r_yz.x
    out$implied$impact  <- imp$impact
    rownames(out$implied) <- NULL
  }
  if (!is.null(out$observed)) {
    out$observed$r_dxj.x  <- sqrt(out$observed$r2dxj.x)
    out$observed$r_yxj.dx <- sqrt(out$observed$r2yxj.dx)
    if (!is.null(out$observed$r2dxj.x_semi))
      out$observed$r_dxj.x_semi <- sqrt(out$observed$r2dxj.x_semi)
    rownames(out$observed) <- NULL
  }
  out
}

#' Impact of the omitted variable on the scale of the ITCV.
#'
#' The ITCV is stated as a product of the omitted variable's partial
#' correlations with the outcome and with the treatment, both given the
#' covariates only. sensemakr's r2yz.dx conditions on the treatment as well,
#' so the outcome correlation is recovered with the recursion formula for
#' partial correlations,
#'   r_yz.x = r_yz.dx * sqrt((1 - r_yd.x^2) (1 - r_dz.x^2)) + r_yd.x * r_dz.x,
#' where r_yd.x is the partial correlation of outcome and treatment given the
#' covariates. sensemakr's bounds identify magnitudes only. The sign of the
#' hypothetical variable is arbitrary, so r_dz.x is taken positive, and
#' r_yz.dx gets the sign of r_yd.x: the omitted variable then biases the
#' estimate away from zero, which is what sensemakr's adjustment with
#' reduce = TRUE removes. With reduce = FALSE the sign is reversed.
#' @noRd
.impact <- function(r2dz.x, r2yz.dx, r_yd.x, reduce = TRUE) {
  r_dz.x  <- sqrt(r2dz.x)
  s       <- if (isTRUE(reduce)) sign(r_yd.x) else -sign(r_yd.x)
  r_yz.dx <- s * sqrt(r2yz.dx)
  r_yz.x  <- r_yz.dx * sqrt((1 - r_yd.x^2) * (1 - r2dz.x)) + r_yd.x * r_dz.x
  list(r_yz.x = r_yz.x, impact = r_yz.x * r_dz.x)
}

#' Attach the implied r2dz.x / r2yz.dx and the impact of each kd row to the
#' results table. The rows are those the adjusted estimates come from: with
#' single benchmarks .get_adjusted_estimates() lines up the manual bound(s)
#' first, then the benchmarks with kd varying within benchmark, and takes
#' the first n_kd rows of that list, or row `bounds_row` with a single kd.
#' Up to 0.4.0.9000 of 2026-11-25 the first benchmark's rows were attached
#' whatever `bounds_row` said and whether or not a manual bound came first.
#' @noRd
.attach_r2 <- function(results, bs, n_kd, r2dz.x = NULL, r2yz.dx = NULL,
                       bounds_row = 1, reduce = TRUE) {
  results$r2dz.x  <- NA_real_
  results$r2yz.dx <- NA_real_
  results$impact  <- NA_real_
  imp <- bs$implied
  if (is.null(imp) || nrow(imp) == 0) return(results)
  if (identical(bs$type, "manual")) {
    results$r2dz.x  <- imp$r2dz.x[1]
    results$r2yz.dx <- imp$r2yz.dx[1]
    results$impact  <- imp$impact[1]
  } else if (identical(bs$type, "single")) {
    all <- imp[, c("r2dz.x", "r2yz.dx", "impact")]
    if (!is.null(r2dz.x)) {
      r2y <- if (is.null(r2yz.dx)) r2dz.x else r2yz.dx
      all <- rbind(data.frame(r2dz.x = r2dz.x, r2yz.dx = r2y,
                              impact = .impact(r2dz.x, r2y, bs$r_yd.x, reduce)$impact),
                   all)
    }
    idx <- if (n_kd > 1) seq_len(n_kd) else as.integer(bounds_row)
    results$r2dz.x  <- all$r2dz.x[idx]
    results$r2yz.dx <- all$r2yz.dx[idx]
    results$impact  <- all$impact[idx]
  } else if (nrow(imp) >= n_kd) {
    results$r2dz.x  <- imp$r2dz.x[seq_len(n_kd)]
    results$r2yz.dx <- imp$r2yz.dx[seq_len(n_kd)]
    results$impact  <- imp$impact[seq_len(n_kd)]
  }
  results
}

#' Print the descriptive benchmark-strength block.
#' @noRd
.print_benchmark_strength <- function(bs) {
  if (is.null(bs)) return(invisible(NULL))
  cat("\nBenchmark strength (descriptive: computed from the data, no standard errors involved)\n")
  if (identical(bs$type, "unavailable")) {
    cat(sprintf("  not available: %s\n", bs$error))
    return(invisible(NULL))
  }
  f4 <- function(v) ifelse(is.na(v), "   .", formatC(v, format = "f", digits = 4))
  f3 <- function(v) ifelse(is.na(v), "   .", formatC(v, format = "f", digits = 3))
  f2 <- function(v) ifelse(is.na(v), "   .", formatC(v, format = "f", digits = 2))

  ob <- bs$observed
  has_semi <- !is.null(ob$r2dxj.x_semi)
  if (!is.null(ob) && nrow(ob) > 0) {
    rows <- unlist(lapply(seq_len(nrow(ob)), function(i) c(
      sprintf("%s with treatment | X", ob$benchmark[i]),
      if (has_semi) sprintf("%s with treatment | X, semi-weights", ob$benchmark[i]),
      sprintf("%s with outcome | D, X", ob$benchmark[i]))))
    r2   <- unlist(lapply(seq_len(nrow(ob)), function(i) c(
      ob$r2dxj.x[i], if (has_semi) ob$r2dxj.x_semi[i], ob$r2yxj.dx[i])))
    tab  <- data.frame(`Partial R2` = f4(r2), `|Partial r|` = f3(sqrt(r2)),
                       check.names = FALSE, stringsAsFactors = FALSE)
    rownames(tab) <- rows
    print(tab, right = TRUE, quote = FALSE)
  }

  imp <- bs$implied
  if (!is.null(imp) && nrow(imp) > 0) {
    cat(if (has_semi)
      "Implied strength of the omitted variable (bounds of Wainstein & Hazlett, 2026, with semi-weights):\n"
      else "Implied strength of the omitted variable (sensemakr bounds):\n")
    tab2 <- data.frame(kd = f2(imp$kd), ky = f2(imp$ky),
                       `R2dz.x` = f4(imp$r2dz.x), `R2yz.dx` = f4(imp$r2yz.dx),
                       Impact = f3(imp$impact),
                       check.names = FALSE, stringsAsFactors = FALSE)
    rownames(tab2) <- if (identical(bs$type, "manual")) imp$benchmark else
      make.unique(sprintf("%sx %s", f2(imp$kd), imp$benchmark), sep = " #")
    print(tab2, right = TRUE, quote = FALSE)
  }
  cat("(R2dz.x, R2yz.dx: partial R2 of the omitted variable with the treatment given the\n",
      "covariates, and with the outcome given treatment and covariates. Impact: product of\n",
      "its partial correlations with the outcome and with the treatment, both given the\n",
      "covariates only, i.e. the scale of the ITCV; signed so that the omitted variable\n",
      "biases the estimate away from zero. Descriptive only: independent of the variance\n",
      "estimator, and not a threshold.)\n", sep = " ")
  invisible(NULL)
}


# ==============================================================================
# Helpers
# ==============================================================================

# Fit the OLS regression, with optional weights. lm() evaluates `weights` in
# `data` and then in the environment of the formula -- not in the function
# that calls lm() -- so a weights vector held in a local variable is not
# found there. The vector is therefore placed in a child of the formula's
# environment, where lm() looks for it.
.fit_lm <- function(formula, data, weight_vec = NULL) {
  if (is.null(weight_vec)) return(lm(formula, data = data))
  env <- new.env(parent = environment(formula))
  assign(".bootmakr_w", weight_vec, envir = env)
  environment(formula) <- env
  lm(formula, data = data, weights = .bootmakr_w)
}
utils::globalVariables(".bootmakr_w")

.resolve_var <- function(x, data, label) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1 && x %in% names(data)) return(data[[x]])
  if (length(x) == nrow(data)) return(x)
  stop(sprintf("`%s` must be a column name in data or vector of length nrow(data).", label))
}

# Parse `cluster` into NULL, list(vec) for one-way, or list(vec1, vec2) for two-way.
.resolve_cluster <- function(cluster, data) {
  if (is.null(cluster)) return(NULL)

  if (is.character(cluster)) {
    if (length(cluster) > 2)
      stop("`cluster` supports at most two dimensions; pass length 1 or 2.")
    bad <- setdiff(cluster, names(data))
    if (length(bad))
      stop(sprintf("Cluster column(s) not found in data: %s",
                   paste(bad, collapse = ", ")))
    return(lapply(cluster, function(nm) data[[nm]]))
  }

  if (is.list(cluster)) {
    if (length(cluster) > 2)
      stop("`cluster` supports at most two dimensions.")
    return(lapply(cluster, function(x) .resolve_var(x, data, "cluster")))
  }

  if (length(cluster) == nrow(data)) return(list(cluster))

  stop("`cluster` must be NULL, column name(s) (length 1 or 2), a list, ",
       "or a vector of length nrow(data).")
}

# Soft warning when a cluster dimension is very small. Cluster-bootstrap
# inference is unreliable with few clusters (Cameron & Miller 2015 rule of
# thumb is 30-50+); at < 10 the CGM two-way variance can degenerate (the
# pathological case is 2: the marginal bootstrap is essentially singular).
.warn_few_clusters <- function(n_clust, label, threshold = 10L) {
  if (!is.null(n_clust) && n_clust < threshold)
    warning(sprintf(
      "Cluster dimension `%s` has only %d unique value(s). The cluster bootstrap on this dimension may give unreliable SEs (rule of thumb: 30+ clusters; <%d is clearly small).",
      label, n_clust, threshold), call. = FALSE)
}

.parse_converge <- function(converge, reps) {
  defaults <- list(do_converge = FALSE, minreps = 500, stepsize = 500,
                   threshold = round(0.75 * reps))
  if (isFALSE(converge)) return(defaults)
  if (isTRUE(converge))  { defaults$do_converge <- TRUE; return(defaults) }
  if (is.list(converge)) {
    opts <- defaults; opts$do_converge <- TRUE
    if (!is.null(converge$minreps))   opts$minreps   <- converge$minreps
    if (!is.null(converge$stepsize))  opts$stepsize   <- converge$stepsize
    if (!is.null(converge$threshold)) opts$threshold  <- converge$threshold
    stopifnot(opts$minreps < reps, opts$stepsize > 0,
              opts$threshold <= reps, opts$threshold >= opts$minreps)
    return(opts)
  }
  stop("`converge` must be TRUE, FALSE, or a list.")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.setup_resampling <- function(data, cluster_vec, strata_vec) {
  n <- nrow(data)
  n_clust <- NULL
  if (!is.null(cluster_vec) && !is.null(strata_vec)) {
    strata_ids <- unique(strata_vec)
    clust_by_strata <- lapply(strata_ids, function(s) {
      idx <- which(strata_vec == s)
      list(cluster_ids = unique(cluster_vec[idx]),
           row_idx = split(idx, cluster_vec[idx]))
    })
    return(list(type = "strat_cluster", n = n,
                n_clust = length(unique(cluster_vec)),
                strata_ids = strata_ids, clust_by_strata = clust_by_strata))
  }
  if (!is.null(cluster_vec)) {
    cluster_ids <- unique(cluster_vec)
    return(list(type = "cluster", n = n, n_clust = length(cluster_ids),
                cluster_ids = cluster_ids, row_idx = split(seq_len(n), cluster_vec)))
  }
  if (!is.null(strata_vec)) {
    return(list(type = "strata", n = n, n_clust = NULL,
                idx_by_strata = split(seq_len(n), strata_vec)))
  }
  list(type = "simple", n = n, n_clust = NULL)
}

.resample_once <- function(info) {
  switch(info$type,
    simple  = sample.int(info$n, replace = TRUE),
    cluster = {
      sampled <- sample(info$cluster_ids, length(info$cluster_ids), replace = TRUE)
      unlist(info$row_idx[as.character(sampled)], use.names = FALSE)
    },
    strata  = {
      idx <- integer(0)
      for (s_idx in info$idx_by_strata) idx <- c(idx, sample(s_idx, length(s_idx), replace = TRUE))
      idx
    },
    strat_cluster = {
      idx <- integer(0)
      for (cs in info$clust_by_strata) {
        sampled <- sample(cs$cluster_ids, length(cs$cluster_ids), replace = TRUE)
        idx <- c(idx, unlist(cs$row_idx[as.character(sampled)], use.names = FALSE))
      }
      idx
    }
  )
}

# Run one bootstrap loop and return the matrix of adjusted estimates plus
# the number of failed reps. Used once for one-way / no-cluster, and three
# times for two-way clustering (G, H, and G x H intersection).
.run_bootstrap_loop <- function(resample_info, reps, data, formula, treat,
                                weight_vec, benchmark_covariates,
                                gbenchmark_covariates, kd, ky, q, alpha,
                                r2dz.x, r2yz.dx, bound_label, reduce,
                                bounds_row, bench_idx = NA_integer_,
                                show_progress, label = NULL,
                                weight_fun = NULL, base_fixed = NULL,
                                semi = NULL) {
  n_kd     <- length(kd)
  boot_mat <- matrix(NA_real_, nrow = reps, ncol = n_kd)
  colnames(boot_mat) <- paste0("kd_", kd)
  # Unadjusted coefficient and its scale in every replicate; see .coef_scale()
  base_mat <- matrix(NA_real_, nrow = reps, ncol = 2,
                     dimnames = list(NULL, c("estimate", "scale")))
  # Partial R2s of the benchmark the table refers to, in every replicate
  # (with semi-weights a third column, the treatment-side partial R2 under
  # the semi-weights)
  bench_cols <- c("r2dxj.x", "r2yxj.dx", if (!is.null(semi)) "r2dxj.x_semi")
  bench_mat <- matrix(NA_real_, nrow = reps, ncol = length(bench_cols),
                      dimnames = list(NULL, bench_cols))
  n_fail   <- 0L

  if (show_progress) {
    clust_msg <- if (!is.null(resample_info$n_clust))
      sprintf(", %d clusters", resample_info$n_clust) else ""
    prefix <- if (!is.null(label)) paste0(label, " ") else ""
    cat(sprintf("%sBootstrapping (%s reps%s)\n", prefix,
                formatC(reps, big.mark = ",", format = "d"), clust_msg))
    pb <- txtProgressBar(min = 0, max = reps, style = 3, width = 50)
  }

  for (b in seq_len(reps)) {
    if (show_progress) setTxtProgressBar(pb, b)
    boot_idx <- .resample_once(resample_info)
    d_boot   <- data[boot_idx, , drop = FALSE]
    # The fit, the unadjusted coefficient and the adjusted estimates are
    # guarded separately: a replicate in which only the benchmark bound fails
    # (kd too large for that resample) still has a coefficient and a scale,
    # and these are all the robustness value needs.
    # Weights: resampled with their units, or re-estimated on the resample
    # when `weights` is a function. A resample on which that function fails
    # is a failed replication.
    w_ok <- TRUE
    w_b  <- if (is.null(weight_fun)) weight_vec[boot_idx] else
      tryCatch(.call_weight_fun(weight_fun, d_boot),
               error = function(e) { w_ok <<- FALSE; NULL })
    fit_b <- if (!w_ok) NULL else
      tryCatch(.fit_lm(formula, d_boot, w_b), error = function(e) NULL)
    if (!is.null(fit_b)) {
      base_mat[b, ] <- tryCatch(.coef_scale(fit_b, treat),
                                error = function(e) c(NA_real_, NA_real_))
      # The benchmark's partial R2s are computed once here and handed on, so
      # that they are kept even when the bound at the requested kd fails.
      # With fix_benchmark they are those of the original fit throughout.
      # Semi-weights: resampled with their units, or re-estimated on the
      # resample when they are functions. Not needed when the benchmark is
      # held fixed. A resample on which a semi-weight function fails keeps
      # its coefficient and scale but has no adjusted estimate.
      semi_ok <- TRUE
      semi_b  <- if (is.null(semi) || !is.null(base_fixed)) NULL else
        tryCatch(.semi_values(semi, d_boot, boot_idx),
                 error = function(e) { semi_ok <<- FALSE; NULL })
      bench_b <- if (!is.null(base_fixed)) base_fixed else if (!semi_ok) NULL else tryCatch(
        .bench_base(fit_b, d_boot, formula, treat, w_b,
                    benchmark_covariates, gbenchmark_covariates, semi_b),
        error = function(e) NULL)
      if (!is.null(bench_b) && !is.na(bench_idx) &&
          bench_idx >= 1L && bench_idx <= length(bench_b$r2dxj.x))
        bench_mat[b, ] <- c(bench_b$r2dxj.x[bench_idx], bench_b$r2yxj.dx[bench_idx],
                            if (!is.null(semi)) bench_b$r2dxj.x_semi[bench_idx])
      if (semi_ok) boot_mat[b, ] <- tryCatch(
        .get_adjusted_estimates(
          fit_b, d_boot, formula, treat, w_b,
          benchmark_covariates, gbenchmark_covariates,
          kd, ky, q, alpha, r2dz.x, r2yz.dx, bound_label, reduce, bounds_row,
          base = bench_b, semi = semi_b
        ), error = function(e) rep(NA_real_, n_kd))
    }
    if (any(!is.finite(boot_mat[b, ]))) n_fail <- n_fail + 1L
  }
  if (show_progress) { close(pb); cat("\n") }

  list(boot_mat = boot_mat, base_mat = base_mat, bench_mat = bench_mat,
       n_fail = n_fail)
}

.compute_boot_stats <- function(boot_mat, obs_estimates, kd, alpha) {
  n_kd   <- length(kd)
  probs  <- c(alpha / 2, 1 - alpha / 2)
  results <- data.frame(kd = kd, estimate = obs_estimates,
                        se = NA_real_, ci_lower = NA_real_,
                        ci_upper = NA_real_, pvalue = NA_real_,
                        stringsAsFactors = FALSE)
  for (i in seq_len(n_kd)) {
    vals <- boot_mat[, i]; vals <- vals[is.finite(vals)]
    if (length(vals) < 10) next
    results$se[i]       <- sd(vals)
    ci                   <- quantile(vals, probs, na.rm = TRUE)
    results$ci_lower[i] <- unname(ci[1])
    results$ci_upper[i] <- unname(ci[2])
    pL <- mean(vals <= 0); pR <- mean(vals >= 0)
    results$pvalue[i] <- 2 * min(pL, pR)
  }
  results
}

# Cameron, Gelbach & Miller (2011) subtractive variance combination for
# two-way cluster bootstrap. Three marginal bootstraps (G, H, intersection)
# yield V_2way = Var_G + Var_H - Var_GH. CIs and p-values use a normal
# approximation since the result is a variance, not a distribution.
# If V_2way < 0 (possible in small samples) we fall back to
# max(Var_G, Var_H) with a warning -- a conservative scalar analogue of the
# eigenvalue truncation used in the multivariate case.
.compute_twoway_stats <- function(boot_samples, obs_estimates, kd, alpha,
                                  crit_df = Inf) {
  n_kd   <- length(kd)
  z_crit <- .crit(alpha, crit_df)
  results <- data.frame(kd = kd, estimate = obs_estimates,
                        se = NA_real_, ci_lower = NA_real_,
                        ci_upper = NA_real_, pvalue = NA_real_,
                        var_G = NA_real_, var_H = NA_real_, var_GH = NA_real_,
                        var_neg_fix = FALSE,
                        stringsAsFactors = FALSE)
  for (i in seq_len(n_kd)) {
    vG  <- boot_samples$G[,  i]; vG  <- vG[is.finite(vG)]
    vH  <- boot_samples$H[,  i]; vH  <- vH[is.finite(vH)]
    vGH <- boot_samples$GH[, i]; vGH <- vGH[is.finite(vGH)]
    if (length(vG) < 10 || length(vH) < 10 || length(vGH) < 10) next

    var_G <- var(vG); var_H <- var(vH); var_GH <- var(vGH)
    var_2 <- var_G + var_H - var_GH

    if (is.finite(var_2) && var_2 < 0) {
      warning(sprintf(
        "Two-way variance estimate negative for kd = %g (V_G + V_H - V_GH = %g). Falling back to max(V_G, V_H) (conservative).",
        kd[i], var_2), call. = FALSE)
      var_2 <- max(var_G, var_H)
      results$var_neg_fix[i] <- TRUE
    }

    se_2 <- sqrt(var_2)
    est  <- obs_estimates[i]
    results$se[i]       <- se_2
    results$ci_lower[i] <- est - z_crit * se_2
    results$ci_upper[i] <- est + z_crit * se_2
    z_stat <- est / se_2
    results$pvalue[i]   <- .pval2(z_stat, crit_df)
    results$var_G[i]    <- var_G
    results$var_H[i]    <- var_H
    results$var_GH[i]   <- var_GH
  }
  results
}

.convergence_diagnostics <- function(boot_vals, obs_estimate, opts) {
  boot_vals <- boot_vals[is.finite(boot_vals)]
  total     <- length(boot_vals)
  reps_seq  <- seq(opts$minreps, total, by = opts$stepsize)
  if (reps_seq[length(reps_seq)] != total) reps_seq <- c(reps_seq, total)
  conv_df <- data.frame(reps = reps_seq, se = NA_real_, pvalue = NA_real_)
  for (j in seq_along(reps_seq)) {
    sub <- boot_vals[seq_len(reps_seq[j])]
    conv_df$se[j] <- sd(sub)
    pL <- mean(sub <= 0); pR <- mean(sub >= 0)
    conv_df$pvalue[j] <- 2 * min(pL, pR)
  }
  thr <- opts$threshold; high <- conv_df$reps >= thr
  se_hi <- if (any(high) && sum(high) > 1) conv_df$se[high] else NA
  p_hi  <- if (any(high) && sum(high) > 1) conv_df$pvalue[high] else NA
  .safe_cv <- function(x) { m <- mean(x); if (m == 0) NA_real_ else sd(x) / m * 100 }
  list(
    data = conv_df,
    summary = list(
      se_mean = mean(conv_df$se), se_range = diff(range(conv_df$se)),
      se_cv = .safe_cv(conv_df$se),
      se_range_hi = if (all(is.na(se_hi))) NA_real_ else diff(range(se_hi)),
      se_cv_hi    = if (all(is.na(se_hi))) NA_real_ else .safe_cv(se_hi),
      p_mean = mean(conv_df$pvalue), p_range = diff(range(conv_df$pvalue)),
      p_cv = .safe_cv(conv_df$pvalue),
      p_range_hi = if (all(is.na(p_hi))) NA_real_ else diff(range(p_hi)),
      p_cv_hi    = if (all(is.na(p_hi))) NA_real_ else .safe_cv(p_hi),
      threshold = thr,
      boot_mean = mean(boot_vals), boot_sd = sd(boot_vals)
    ),
    boot_vals = boot_vals, obs_estimate = obs_estimate, opts = opts
  )
}

# Convergence diagnostics for the CGM two-way bootstrap.
# Tracks the combined SE and normal-approximation p-value at each cumulative
# rep count using the same V_G + V_H - V_GH formula as the headline result.
# The intersection bootstrap distribution is exposed as boot_vals so the
# existing convergence plot can draw a histogram.
.convergence_diagnostics_twoway <- function(boot_samples, obs_estimate, alpha, opts,
                                            crit_df = Inf) {
  vG  <- boot_samples$G[,  1]
  vH  <- boot_samples$H[,  1]
  vGH <- boot_samples$GH[, 1]
  total    <- min(sum(is.finite(vG)), sum(is.finite(vH)), sum(is.finite(vGH)))
  reps_seq <- seq(opts$minreps, total, by = opts$stepsize)
  if (reps_seq[length(reps_seq)] != total) reps_seq <- c(reps_seq, total)
  conv_df  <- data.frame(reps = reps_seq, se = NA_real_, pvalue = NA_real_)
  for (j in seq_along(reps_seq)) {
    n <- reps_seq[j]
    v_G  <- var(vG[seq_len(n)],  na.rm = TRUE)
    v_H  <- var(vH[seq_len(n)],  na.rm = TRUE)
    v_GH <- var(vGH[seq_len(n)], na.rm = TRUE)
    v2   <- v_G + v_H - v_GH
    if (!is.finite(v2) || v2 < 0) v2 <- max(v_G, v_H, na.rm = TRUE)
    se   <- sqrt(v2)
    conv_df$se[j]     <- se
    conv_df$pvalue[j] <- .pval2(obs_estimate / se, crit_df)
  }
  thr <- opts$threshold; high <- conv_df$reps >= thr
  se_hi <- if (any(high) && sum(high) > 1) conv_df$se[high] else NA
  p_hi  <- if (any(high) && sum(high) > 1) conv_df$pvalue[high] else NA
  .safe_cv <- function(x) { m <- mean(x); if (m == 0) NA_real_ else sd(x) / m * 100 }
  vGH_clean <- vGH[is.finite(vGH)]
  list(
    data = conv_df,
    summary = list(
      se_mean = mean(conv_df$se), se_range = diff(range(conv_df$se)),
      se_cv = .safe_cv(conv_df$se),
      se_range_hi = if (all(is.na(se_hi))) NA_real_ else diff(range(se_hi)),
      se_cv_hi    = if (all(is.na(se_hi))) NA_real_ else .safe_cv(se_hi),
      p_mean = mean(conv_df$pvalue), p_range = diff(range(conv_df$pvalue)),
      p_cv = .safe_cv(conv_df$pvalue),
      p_range_hi = if (all(is.na(p_hi))) NA_real_ else diff(range(p_hi)),
      p_cv_hi    = if (all(is.na(p_hi))) NA_real_ else .safe_cv(p_hi),
      threshold = thr,
      boot_mean = mean(vGH_clean), boot_sd = sd(vGH_clean)
    ),
    boot_vals = vGH_clean, obs_estimate = obs_estimate, opts = opts,
    twoway = TRUE
  )
}


# ==============================================================================
# Breakdown kd
# ==============================================================================
# The smallest kd (with ky = kd) at which the bootstrap test of the adjusted
# estimate stops rejecting. sensemakr's bound on the omitted variable is a
# closed form in kd and the benchmark's two partial R2s, and the adjusted
# estimate a closed form in that bound, the coefficient and its scale. With
# these four numbers stored for every replicate, the adjusted estimate of
# every replicate can be evaluated at any kd without refitting. With
# semi-weights there is a fifth number, the benchmark's partial R2 with the
# treatment under the semi-weights (`r2s` below), and the bound is that of
# Wainstein and Hazlett (2026); with r2s equal to r2d it is sensemakr's.

#' Partial R2s of the benchmark(s) for one fit: the group if
#' gbenchmark_covariates is given, otherwise one entry per benchmark
#' covariate; NULL without a benchmark.
#' @noRd
.bench_base <- function(fit, data, formula, treat, weight_vec,
                        benchmark_covariates, gbenchmark_covariates, semi = NULL) {
  if (!is.null(gbenchmark_covariates))
    return(.group_base_r2(fit, data, formula, treat, weight_vec, gbenchmark_covariates,
                          semi))
  if (!is.null(benchmark_covariates))
    return(.single_base_r2(fit, treat, benchmark_covariates, semi))
  NULL
}

#' Index of the benchmark that the table of adjusted estimates refers to.
#' The rows of sensemakr's bounds table are: manual bound(s) first, then the
#' benchmarks, kd varying within benchmark; the table takes the first n_kd
#' rows of it, or row `bounds_row` when there is a single kd. NA when the
#' table reports a manual row or there is no benchmark.
#' @noRd
.table_benchmark_index <- function(benchmark_covariates, gbenchmark_covariates,
                                   r2dz.x, n_kd, bounds_row) {
  if (!is.null(gbenchmark_covariates)) return(1L)
  if (is.null(benchmark_covariates))   return(NA_integer_)
  n_man <- length(r2dz.x)
  if (n_kd > 1) return(if (n_man == 0) 1L else NA_integer_)
  j <- as.integer(bounds_row) - n_man
  if (is.na(j) || j < 1L || j > length(benchmark_covariates)) NA_integer_ else j
}

#' Largest kd the benchmark allows: sensemakr's bound on r2dz.x reaches one,
#' and the bound stops with an error, at kd = (1 - r2dxj.x) / r2dxj.x. With
#' semi-weights the bound on r2dz.x is kd * r2s / (1 - r2d), so the limit is
#' (1 - r2d) / r2s; the second condition of the bound, r2zxj.xd < 1, is
#' kd * r2s < 1 - r2d as well.
#' @noRd
.kd_limit <- function(r2d, r2s = r2d) ifelse(r2s > 0, (1 - r2d) / r2s, Inf)

#' Bound on the two partial R2s of the omitted variable, vectorised over
#' benchmarks (replicates) and/or kd. `r2d` and `r2y` are the benchmark's
#' partial R2 with the treatment and with the outcome under the weights of
#' the regression, `r2s` its partial R2 with the treatment under the
#' semi-weights. Wainstein and Hazlett (2026, equations 25 and 26):
#'   r2dz.x  = kd * r2s / (1 - r2d)
#'   r2zxj   = kd * r2s / (1 - kd * r2s) * r2d / (1 - r2d)
#'   r2yz.dx = ((sqrt(ky) + sqrt(r2zxj)) / sqrt(1 - r2zxj))^2 * r2y / (1 - r2y)
#' With r2s = r2d (the default: no semi-weights) this is sensemakr's
#' ovb_partial_r2_bound(). NA where sensemakr stops with an error (kd at or
#' beyond the limit); r2yz.dx is capped at one, as sensemakr does, and
#' `capped` says where.
#' @noRd
.bound_at <- function(r2d, r2y, k, r2s = r2d, ky = k) {
  force(r2s); force(ky)
  ok   <- is.finite(r2d) & is.finite(r2y) & is.finite(r2s) & r2d >= 0 & r2d < 1 &
    r2y >= 0 & r2y <= 1 & r2s >= 0 & r2s <= 1 & k >= 0 & k < .kd_limit(r2d, r2s)
  r2d  <- ifelse(ok, r2d, NA_real_)
  r2dz <- k * r2s / (1 - r2d)
  r2zx <- k * (r2s * r2d) / ((1 - k * r2s) * (1 - r2d))
  bad  <- !is.na(r2dz) & (r2dz >= 1 | r2zx >= 1 | r2zx < 0)
  r2dz[bad] <- NA_real_; r2zx[bad] <- NA_real_
  raw  <- ((sqrt(ky) + sqrt(r2zx)) / sqrt(1 - r2zx))^2 * (r2y / (1 - r2y))
  list(r2dz.x = r2dz, r2yz.dx = pmin(raw, 1), capped = !is.na(raw) & raw > 1)
}

#' Bound for one benchmark, as a row per kd. Without semi-weights this is
#' sensemakr::ovb_partial_r2_bound() itself, so that results without
#' `semiweights` are those of sensemakr to the last digit. With them
#' (`r2dxj.x_semi`, the benchmark's partial R2 with the treatment under the
#' semi-weights) it is .bound_at(), with sensemakr's error when kd is at or
#' beyond the limit and its warning when r2yz.dx is capped at one.
#' @noRd
.ovb_bound <- function(r2dxj.x, r2yxj.dx, kd = 1, ky = kd, bound_label = "manual",
                       r2dxj.x_semi = NULL) {
  if (is.null(r2dxj.x_semi))
    return(sensemakr::ovb_partial_r2_bound(
      r2dxj.x = r2dxj.x, r2yxj.dx = r2yxj.dx, kd = kd, ky = ky,
      bound_label = bound_label))
  r2 <- c(r2dxj.x, r2yxj.dx, r2dxj.x_semi)
  if (!all(is.finite(r2)) || any(r2 < 0 | r2 > 1))
    stop("partial R2 must be a number between zero and one")
  b <- .bound_at(r2dxj.x, r2yxj.dx, kd, r2dxj.x_semi, ky)
  if (anyNA(b$r2dz.x))
    stop("Implied bound on r2dz.x >= 1. Impossible kd value. Try a lower kd.")
  if (any(b$capped))
    warning("Implied bound on r2yz.dx greater than 1, try lower kd and/or ky. Setting r2yz.dx to 1.")
  data.frame(bound_label = bound_label, r2dz.x = b$r2dz.x, r2yz.dx = b$r2yz.dx,
             stringsAsFactors = FALSE)
}

#' Partial R2 of the benchmark with the treatment under the semi-weights in
#' a matrix of stored benchmark R2s: its third column if there is one,
#' otherwise the first (no semi-weights).
#' @noRd
.bench_semi <- function(bench) if (ncol(bench) >= 3) bench[, 3] else bench[, 1]

#' Adjusted estimate (reduce = TRUE) for an omitted variable k times as
#' strong as the benchmark; NA where the bound does not exist.
#' @noRd
.adj_kd <- function(est, scale, r2d, r2y, k, r2s = r2d) {
  b <- .bound_at(r2d, r2y, k, r2s)
  sign(est) * (abs(est) - scale * sqrt(b$r2yz.dx * b$r2dz.x / (1 - b$r2dz.x)))
}

#' kd at which each adjusted estimate reaches zero, by bisection. The bias
#' increases in kd and is unbounded at the limit, so there is one crossing
#' below the limit; the upper end of the final bracket is returned, where the
#' adjusted estimate is at or just beyond zero. Inf if the bias is zero
#' throughout (benchmark unrelated to treatment or outcome).
#' @noRd
.kd_zero <- function(est, scale, r2d, r2y, r2s = r2d) {
  lim <- .kd_limit(r2d, r2s)
  lo  <- rep(0, length(est))
  hi  <- ifelse(is.finite(lim), lim * (1 - 1e-12), NA_real_)
  g   <- function(k) abs(est) - scale * sqrt({
    b <- .bound_at(r2d, r2y, k, r2s); b$r2yz.dx * b$r2dz.x / (1 - b$r2dz.x) })
  never <- is.na(hi) | !is.finite(g(hi)) | g(hi) > 0
  hi[never] <- 0
  for (i in seq_len(200)) {
    mid  <- (lo + hi) / 2
    done <- mid <= lo | mid >= hi
    if (all(done | never)) break
    neg  <- g(mid) <= 0
    hi   <- ifelse(!done & neg, mid, hi)
    lo   <- ifelse(!done & !neg, mid, lo)
  }
  hi[never] <- Inf
  hi
}

#' Breakdown kd, percentile case.
#'
#' The smallest kd from which the percentile p-value of the table,
#' 2 * min(share <= 0, share >= 0) over the replicates that have a bound at
#' that kd, is at least alpha. The p-value changes only where a replicate's
#' adjusted estimate crosses zero or where a replicate reaches its own limit
#' and drops out, so the search over those points is exact. As in the table,
#' the adjustment follows the sign of the coefficient in each replicate, and
#' fewer than 10 usable replicates give no p-value. The search ends at
#' `k_limit`, the limit of the original fit.
#' @noRd
.breakdown_percentile <- function(base, bench, alpha, k_limit) {
  semi_col <- .bench_semi(bench)
  ok <- is.finite(base[, "estimate"]) & is.finite(base[, "scale"]) & base[, "scale"] > 0 &
    is.finite(bench[, 1]) & is.finite(bench[, 2]) &
    bench[, 1] >= 0 & bench[, 1] < 1 & bench[, 2] >= 0 & bench[, 2] <= 1 &
    is.finite(semi_col) & semi_col >= 0 & semi_col <= 1
  e <- base[ok, "estimate"]; s <- base[ok, "scale"]
  r2d <- bench[ok, 1]; r2y <- bench[ok, 2]; r2s <- semi_col[ok]; B <- length(e)
  out <- list(kd = NA_real_, n_used = NA_integer_, n_reps = B, note = NULL)
  if (B < 10) { out$note <- "Fewer than 10 usable replicates."; return(out) }
  if (2 * min(mean(e <= 0), mean(e >= 0)) >= alpha) { out$kd <- 0; out$n_used <- B; return(out) }

  lim   <- .kd_limit(r2d, r2s)
  cross <- pmin(.kd_zero(e, s, r2d, r2y, r2s), lim) # crossing, or dropping out
  cand  <- sort(unique(c(cross, lim)))
  cand  <- cand[is.finite(cand) & cand < k_limit]
  if (length(cand) == 0) {
    out$note <- "The estimate stays significant for every kd the benchmark allows."
    return(out)
  }
  cnt <- function(x) findInterval(cand, sort(x))    # number of x at or below cand
  pos <- e > 0; neg <- e < 0; zer <- e == 0
  n_valid <- B - cnt(lim)
  n_le0 <- (cnt(cross[pos]) - cnt(lim[pos])) + (sum(neg) - cnt(cross[neg])) + (sum(zer) - cnt(lim[zer]))
  n_ge0 <- (sum(pos) - cnt(cross[pos])) + (cnt(cross[neg]) - cnt(lim[neg])) + (sum(zer) - cnt(lim[zer]))
  pval  <- ifelse(n_valid >= 10, 2 * pmin(n_le0 / n_valid, n_ge0 / n_valid), NA_real_)
  i <- which(is.na(pval) | pval >= alpha)[1]
  if (is.na(i)) {
    out$note <- "The estimate stays significant for every kd the benchmark allows."
  } else if (is.na(pval[i])) {
    out$note <- "Too few replicates have a bound at the kd where the estimate might stop being significant."
  } else {
    out$kd <- cand[i]; out$n_used <- as.integer(n_valid[i])
  }
  out
}

#' First point on [0, upper] at which h() is at or below zero, for a
#' function that starts above zero: located on a grid, then refined with
#' uniroot inside the first bracket. NA if h() is missing before that point
#' or never gets there.
#' @noRd
.first_crossing <- function(h, upper, n = 401) {
  grid <- seq(0, upper, length.out = n)
  hv   <- vapply(grid, h, numeric(1))
  i    <- which(is.na(hv) | hv <= 0)[1]
  if (is.na(i) || is.na(hv[i]) || i == 1) return(if (!is.na(i) && i == 1 && !is.na(hv[1])) 0 else NA_real_)
  tryCatch(stats::uniroot(h, lower = grid[i - 1], upper = grid[i], tol = 1e-12)$root,
           error = function(e) NA_real_)
}

#' Breakdown kd, two-way case: the smallest kd at which
#' |adjusted estimate| / sqrt(V_G + V_H - V_GH) falls to the normal critical
#' value, the variances taken over the replicates that have a bound at that
#' kd (at least 10 in each of the three bootstraps).
#' @noRd
.breakdown_twoway <- function(base_list, bench_list, obs, obs_bench, alpha,
                              crit_df = Inf) {
  out <- list(kd = NA_real_, n_used = NA_integer_, n_reps = NA_integer_, note = NULL)
  dims <- c("G", "H", "GH")
  adj  <- function(d, k) .adj_kd(base_list[[d]][, "estimate"], base_list[[d]][, "scale"],
                                 bench_list[[d]][, 1], bench_list[[d]][, 2], k,
                                 .bench_semi(bench_list[[d]]))
  out$n_reps <- min(vapply(dims, function(d) sum(is.finite(adj(d, 0))), integer(1)))
  z_crit <- .crit(alpha, crit_df)
  h <- function(k) {
    v <- vapply(dims, function(d) {
      a <- adj(d, k); a <- a[is.finite(a)]
      if (length(a) < 10) NA_real_ else var(a)
    }, numeric(1))
    if (anyNA(v)) return(NA_real_)
    v2 <- v[["G"]] + v[["H"]] - v[["GH"]]
    if (!is.finite(v2) || v2 < 0) v2 <- max(v[["G"]], v[["H"]])
    abs(.adj_kd(obs[[1]], obs[[2]], obs_bench[1], obs_bench[2], k, obs_bench[3])) /
      sqrt(v2) - z_crit
  }
  # the adjusted estimate of the original fit is zero at its own crossing
  upper <- .kd_zero(obs[[1]], obs[[2]], obs_bench[1], obs_bench[2], obs_bench[3])
  if (!is.finite(upper)) {
    out$note <- "The estimate stays significant for every kd the benchmark allows."
    return(out)
  }
  out$kd <- .first_crossing(h, upper)
  if (is.na(out$kd)) {
    out$note <- "Too few replicates have a bound at the kd where the estimate might stop being significant."
  } else {
    out$n_used <- min(vapply(dims, function(d) sum(is.finite(adj(d, out$kd))), integer(1)))
  }
  out
}

#' Breakdown kd under the conventional standard error: the smallest kd at
#' which the confidence interval that sensemakr::ovb_bounds() reports for
#' the adjusted estimate (adjusted_se, t critical value with the residual
#' degrees of freedom) includes zero. The adjusted standard error shrinks as
#' r2yz.dx grows, so the t-statistic need not be monotone in kd; the first
#' crossing is taken.
#' @noRd
.breakdown_conventional <- function(obs, obs_bench, dof, alpha) {
  upper <- .kd_zero(obs[[1]], obs[[2]], obs_bench[1], obs_bench[2], obs_bench[3])
  if (!is.finite(upper) || !is.finite(dof) || dof < 2) return(NA_real_)
  t_crit <- stats::qt(1 - alpha / 2, df = dof)
  h <- function(k) {
    b  <- .bound_at(obs_bench[1], obs_bench[2], k, obs_bench[3])
    se <- sqrt((1 - b$r2yz.dx) / (1 - b$r2dz.x)) * obs[[2]] / sqrt(dof) * sqrt(dof / (dof - 1))
    a  <- .adj_kd(obs[[1]], obs[[2]], obs_bench[1], obs_bench[2], k, obs_bench[3])
    if (!is.finite(se) || se <= 0) return(if (isTRUE(a == 0)) -t_crit else Inf)
    abs(a) / se - t_crit
  }
  .first_crossing(h, upper)
}

#' Breakdown kd of the benchmark the table refers to
#' @noRd
.breakdown <- function(fit, treat, alpha, reduce, kd, ky, bench_idx,
                       bench_strength, boot_base, boot_bench, is_twoway,
                       crit_df = Inf) {
  out <- list(kd = NA_real_, kd_conventional = NA_real_, kd_limit = NA_real_,
              benchmark = NA_character_, alpha = alpha,
              method = if (is_twoway) "cgm_twoway" else "percentile",
              n_used = NA_integer_, n_reps = NA_integer_, note = NULL)
  obs_tab <- bench_strength$observed
  if (is.na(bench_idx) || is.null(obs_tab) || bench_idx > nrow(obs_tab)) {
    out$note <- "A breakdown kd needs a benchmark, and the table must report it rather than a manual R2."
    return(out)
  }
  out$benchmark <- obs_tab$benchmark[bench_idx]
  # third entry: the partial R2 with the treatment under the semi-weights,
  # or under the weights again when there are none
  obs_bench <- c(obs_tab$r2dxj.x[bench_idx], obs_tab$r2yxj.dx[bench_idx],
                 if (is.null(obs_tab$r2dxj.x_semi)) obs_tab$r2dxj.x[bench_idx]
                 else obs_tab$r2dxj.x_semi[bench_idx])
  out$kd_limit <- unname(.kd_limit(obs_bench[1], obs_bench[3]))
  if (!isTRUE(reduce)) {
    out$note <- "The breakdown kd is defined for reduce = TRUE only."
    return(out)
  }
  obs <- tryCatch(.coef_scale(fit, treat),
                  error = function(e) c(estimate = NA_real_, scale = NA_real_))
  if (!all(is.finite(obs)) || !all(is.finite(obs_bench))) {
    out$note <- "The breakdown kd could not be computed."
    return(out)
  }
  out$kd_conventional <- tryCatch(
    .breakdown_conventional(obs, obs_bench, fit$df.residual, alpha),
    error = function(e) NA_real_)
  res <- tryCatch(
    if (is_twoway) .breakdown_twoway(boot_base, boot_bench, obs, obs_bench, alpha, crit_df)
    else           .breakdown_percentile(boot_base, boot_bench, alpha, out$kd_limit),
    error = function(e) list(kd = NA_real_, n_used = NA_integer_, n_reps = NA_integer_,
                             note = "The breakdown kd could not be computed."))
  out$kd <- res$kd; out$n_used <- res$n_used; out$n_reps <- res$n_reps
  notes <- res$note
  if (!isTRUE(all.equal(as.numeric(kd), as.numeric(ky))))
    notes <- c(notes, "The breakdown kd takes ky = kd; the call uses a different ky.")
  if (!is.na(out$kd) && !is.na(out$n_used) && out$n_used < out$n_reps)
    notes <- c(notes, sprintf(
      "At this kd %d of %d replicates have no bound and are left out, as in the table.",
      out$n_reps - out$n_used, out$n_reps))
  out$note <- if (length(notes)) paste(notes, collapse = " ") else NULL
  out
}

.print_breakdown <- function(bd) {
  if (is.null(bd) || is.na(bd$benchmark)) return(invisible(NULL))
  pct <- format(bd$alpha * 100)
  f   <- function(v) if (is.na(v)) "    NA" else formatC(v, format = "f", digits = 4)
  cat(sprintf("\nBreakdown kd (smallest kd at which the estimate is not significant at %s%%;\n", pct),
      sprintf(" benchmark: %s, ky = kd):\n", bd$benchmark), sep = "")
  lab <- c("Conventional SE", "Bootstrap")
  cat(sprintf("  %-*s  %s\n", max(nchar(lab)), lab, c(f(bd$kd_conventional), f(bd$kd))), sep = "")
  if (is.finite(bd$kd_limit))
    cat(sprintf("(largest kd this benchmark allows: %s)\n", formatC(bd$kd_limit, format = "f", digits = 2)))
  if (!is.null(bd$note)) cat(strwrap(paste("Note:", bd$note), width = 78), sep = "\n")
  invisible(NULL)
}


# ==============================================================================
# Robustness value
# ==============================================================================
# With an omitted variable of strength r2dz.x = r2yz.dx = r, sensemakr's
# adjusted estimate is
#     sign(est) * (|est| - scale * r / sqrt(1 - r)),   scale = se * sqrt(dof),
# so the coefficient and its scale are all that is needed to evaluate it at
# any r, in the original fit and in every bootstrap replicate, without
# refitting anything. `scale` is the ratio of the residual standard
# deviations of the outcome and of the treatment; the conventional standard
# error enters it only as a way of computing that ratio.

#' Coefficient of the treatment and its scale, se * sqrt(dof)
#' @noRd
.coef_scale <- function(fit, treat) {
  coefs <- summary(fit)$coefficients
  if (!(treat %in% rownames(coefs)))
    return(c(estimate = NA_real_, scale = NA_real_))
  c(estimate = unname(coefs[treat, "Estimate"]),
    scale    = unname(coefs[treat, "Std. Error"]) * sqrt(fit$df.residual))
}

#' sensemakr's adjusted estimate (reduce = TRUE) at equal strength r
#' @noRd
.adj_at <- function(est, scale, r) sign(est) * (abs(est) - scale * r / sqrt(1 - r))

#' Strength r at which the adjusted estimate is exactly zero. This is
#' sensemakr's RV_q for q = 1, written so that it does not lose precision
#' when |est| / scale is large.
#' @noRd
.rv_zero <- function(est, scale) {
  f2 <- (est / scale)^2
  2 * f2 / (sqrt(f2^2 + 4 * f2) + f2)
}

#' Bootstrap robustness value, percentile case.
#'
#' The smallest r from which the percentile-bootstrap p-value of the adjusted
#' estimate, 2 * min(share <= 0, share >= 0), is at least alpha. The p-value
#' is a step function of r that changes only where a replicate's adjusted
#' estimate crosses zero, i.e. at that replicate's own .rv_zero(), so the
#' search over those points is exact. As in the main table, the direction of
#' the adjustment follows the sign of the coefficient in each replicate: a
#' replicate with a positive coefficient is at or below zero from its
#' crossing point onwards, one with a negative coefficient up to it. Both
#' groups are counted, because the p-value need not be monotone in r.
#' @noRd
.bootstrap_rv_percentile <- function(base, alpha) {
  ok <- is.finite(base[, 1]) & is.finite(base[, 2]) & base[, 2] > 0
  e  <- base[ok, 1]; s <- base[ok, 2]; B <- length(e)
  if (B < 10) return(NA_real_)
  if (2 * min(mean(e <= 0), mean(e >= 0)) >= alpha) return(0)

  cr   <- .rv_zero(e, s)
  cpos <- sort(cr[e > 0]); cneg <- sort(cr[e < 0]); n0 <- sum(e == 0)
  cand <- sort(unique(cr[e != 0]))
  pos_le <- findInterval(cand, cpos)      # positive replicates crossed by cand
  neg_le <- findInterval(cand, cneg)      # negative replicates crossed by cand
  # counts just above each candidate point
  n_le0 <- pos_le + (length(cneg) - neg_le) + n0
  n_ge0 <- (length(cpos) - pos_le) + neg_le + n0
  hit   <- which(2 * pmin(n_le0 / B, n_ge0 / B) >= alpha)
  if (length(hit) == 0) NA_real_ else cand[hit[1]]
}

#' Bootstrap robustness value, two-way case.
#'
#' The smallest r at which |adjusted estimate| / SE falls to the normal
#' critical value, with SE = sqrt(V_G + V_H - V_GH) recomputed at every r
#' from the three sets of replicates (and the same max(V_G, V_H) fallback as
#' the main table, without its warning). A grid locates the first crossing;
#' uniroot then refines it inside that bracket.
#' @noRd
.bootstrap_rv_twoway <- function(base_list, obs, alpha, crit_df = Inf) {
  bs <- lapply(base_list, function(m)
    m[is.finite(m[, 1]) & is.finite(m[, 2]) & m[, 2] > 0, , drop = FALSE])
  if (any(vapply(bs, nrow, integer(1)) < 10) || !all(is.finite(obs)))
    return(NA_real_)
  z_crit <- .crit(alpha, crit_df)
  h <- function(r) {
    v  <- vapply(bs, function(m) var(.adj_at(m[, 1], m[, 2], r)), numeric(1))
    v2 <- v[["G"]] + v[["H"]] - v[["GH"]]
    if (!is.finite(v2) || v2 < 0) v2 <- max(v[["G"]], v[["H"]])
    abs(.adj_at(obs[[1]], obs[[2]], r)) / sqrt(v2) - z_crit
  }
  if (h(0) <= 0) return(0)
  # the adjusted estimate of the original fit is zero at .rv_zero(), so the
  # statistic has fallen below the critical value somewhere before it
  grid <- seq(0, .rv_zero(obs[[1]], obs[[2]]), length.out = 201)
  hv   <- vapply(grid, h, numeric(1))
  i    <- which(hv <= 0)[1]
  if (is.na(i)) return(NA_real_)
  stats::uniroot(h, lower = grid[i - 1], upper = grid[i], tol = 1e-12)$root
}

#' Robustness values of the treatment coefficient
#'
#' rv_q and rv_qa are sensemakr's own numbers for the original fit; rv_qa
#' rests on the conventional OLS standard error. rv_qa_boot is the bootstrap
#' counterpart: the strength at which the test in the main table stops
#' rejecting. It is defined for q = 1 and reduce = TRUE only, because that
#' is the only null hypothesis the bootstrap table tests (adjusted
#' estimate = 0, omitted variable pulling the estimate towards zero).
#' @noRd
.robustness <- function(fit, treat, q, alpha, reduce, boot_base, is_twoway,
                        crit_df = Inf) {
  sm_rv <- function(a) tryCatch(
    as.numeric(sensemakr::robustness_value(fit, covariates = treat, q = q, alpha = a)),
    error = function(e) NA_real_)
  obs  <- tryCatch(.coef_scale(fit, treat),
                   error = function(e) c(estimate = NA_real_, scale = NA_real_))
  note <- NULL
  rv_boot <- NA_real_
  if (!isTRUE(reduce)) {
    note <- "The bootstrap robustness value is defined for reduce = TRUE only."
  } else if (!isTRUE(all.equal(q, 1))) {
    note <- "The bootstrap robustness value is available for q = 1 only."
  } else {
    rv_boot <- tryCatch(
      if (is_twoway) .bootstrap_rv_twoway(boot_base, obs, alpha, crit_df)
      else           .bootstrap_rv_percentile(boot_base, alpha),
      error = function(e) NA_real_)
    if (is.na(rv_boot))
      note <- "The bootstrap robustness value could not be computed."
  }
  n_used <- if (is_twoway) {
    min(vapply(boot_base, function(m) sum(is.finite(m[, 1]) & is.finite(m[, 2])), integer(1)))
  } else sum(is.finite(boot_base[, 1]) & is.finite(boot_base[, 2]))
  list(rv_q = sm_rv(1), rv_qa = sm_rv(alpha), rv_qa_boot = rv_boot,
       q = q, alpha = alpha,
       method = if (is_twoway) "cgm_twoway" else "percentile",
       n_used = n_used, estimate = obs[["estimate"]], scale = obs[["scale"]],
       note = note)
}

.print_robustness <- function(rb, crit_df = Inf) {
  if (is.null(rb)) return(invisible(NULL))
  pct <- format(rb$alpha * 100)
  f   <- function(v) if (is.na(v)) "    NA" else formatC(v, format = "f", digits = 4)
  cat("\nRobustness value (omitted variable equally strong with treatment and outcome):\n")
  lab <- c(sprintf("Adjusted estimate = %s", if (isTRUE(all.equal(rb$q, 1))) "0"
                   else sprintf("%g%% lower", rb$q * 100)),
           sprintf("Not significant at %s%%, conventional SE (sensemakr)", pct),
           sprintf("Not significant at %s%%, bootstrap", pct))
  val <- c(f(rb$rv_q), f(rb$rv_qa), f(rb$rv_qa_boot))
  cat(sprintf("  %-*s  %s\n", max(nchar(lab)), lab, val), sep = "")
  cat("(partial R2 of the omitted variable with the treatment and with the outcome;\n",
      if (identical(rb$method, "cgm_twoway") && is.finite(crit_df))
        " bootstrap: smallest strength at which the two-way t-test stops rejecting)\n"
      else if (identical(rb$method, "cgm_twoway"))
        " bootstrap: smallest strength at which the two-way z-test stops rejecting)\n"
      else
        " bootstrap: smallest strength at which the bootstrap p-value reaches alpha)\n",
      sep = "")
  if (!is.null(rb$note)) cat("Note: ", rb$note, "\n", sep = "")
  invisible(NULL)
}


# ==============================================================================
# Critical values and weight functions
# ==============================================================================

#' Two-sided critical value and p-value of the two-way test: normal, or t
#' when `df` is finite.
#' @noRd
.crit  <- function(alpha, df = Inf)
  if (is.finite(df)) stats::qt(1 - alpha / 2, df = df) else qnorm(1 - alpha / 2)
.pval2 <- function(stat, df = Inf)
  if (is.finite(df)) 2 * stats::pt(-abs(stat), df = df) else 2 * pnorm(-abs(stat))

#' Call a user-supplied weight function on a data set and check what it
#' returns.
#' @noRd
.call_weight_fun <- function(fun, data, label = "weights") {
  w <- fun(data)
  if (!is.numeric(w) || length(w) != nrow(data) || anyNA(w) || any(w < 0))
    stop("A `", label, "` function must return one non-negative number per row of ",
         "the data it is given.", call. = FALSE)
  as.numeric(w)
}

#' Check `semiweights` and bring it into one form: NULL, or a list with
#' `kind` ("vector" or "function") and `items`, one per benchmark covariate
#' (one for a grouped benchmark). Semi-weights must be of the same kind as
#' the weights, so that both are resampled with their units or both are
#' re-estimated on every resample.
#' @noRd
.resolve_semiweights <- function(semiweights, data, weights, weight_fun,
                                 benchmark_covariates, gbenchmark_covariates) {
  if (is.null(semiweights)) return(NULL)
  if (is.null(weights))
    stop("`semiweights` needs `weights`: semi-weights are the weights formed ",
         "without the benchmark.", call. = FALSE)
  n_b <- if (!is.null(gbenchmark_covariates)) 1L else length(benchmark_covariates)
  if (n_b == 0L)
    stop("`semiweights` needs a benchmark (`benchmark_covariates` or ",
         "`gbenchmark_covariates`).", call. = FALSE)
  items <- if (is.list(semiweights) && !is.data.frame(semiweights)) semiweights
           else if (is.character(semiweights) && all(semiweights %in% names(data)))
             as.list(semiweights)
           else list(semiweights)
  if (length(items) != n_b)
    stop(sprintf(paste0(
      "`semiweights` has %d entr%s but there %s %d benchmark%s: give a list with one ",
      "entry per benchmark covariate (a grouped benchmark takes one)."),
      length(items), if (length(items) == 1) "y" else "ies",
      if (n_b == 1) "is" else "are", n_b, if (n_b == 1) "" else "s"), call. = FALSE)
  is_fun <- vapply(items, is.function, logical(1))
  if (!is.null(weight_fun)) {
    if (!all(is_fun))
      stop("`weights` is a function, so `semiweights` must be a function too (or a ",
           "list of functions): both are then re-estimated on every resample.",
           call. = FALSE)
    return(list(kind = "function", items = items))
  }
  if (any(is_fun))
    stop("`weights` is a column or a vector, so `semiweights` must be one too (or a ",
         "list of them): both are then resampled with their units.", call. = FALSE)
  vecs <- lapply(items, function(v) {
    v <- .resolve_var(v, data, "semiweights")
    if (!is.numeric(v) || anyNA(v) || any(v < 0))
      stop("`semiweights` must be non-negative numbers without missing values.",
           call. = FALSE)
    as.numeric(v)
  })
  list(kind = "vector", items = vecs)
}

#' Semi-weights for one data set, a list with one vector per benchmark:
#' functions are called on `data`; vectors are returned as they are, or the
#' entries `idx` of them for a resample. NULL without semi-weights.
#' @noRd
.semi_values <- function(semi, data, idx = NULL) {
  if (is.null(semi)) return(NULL)
  if (identical(semi$kind, "function"))
    return(lapply(semi$items, .call_weight_fun, data = data, label = "semiweights"))
  if (is.null(idx)) semi$items else lapply(semi$items, function(v) v[idx])
}


# ==============================================================================
# Standard errors side by side, and the robust-to-conventional ratio
# ==============================================================================
# The table of adjusted estimates rests on the bootstrap. Two closed-form
# standard errors exist for the same adjusted estimate and are reported next
# to it: sensemakr's, built on the conventional OLS standard error, and the
# rescaled robust standard error of Wainstein and Hazlett (2026, Option 3 and
# Appendix B.2). Both refer to the regression that would include the omitted
# variable; the bootstrap refers to the sampling variability of the adjusted
# estimate. The robust standard errors are computed here from the scores of
# the treatment coefficient (Frisch-Waugh-Lovell), so that no further
# package is needed: HC1 without clusters, CR1 with one cluster variable,
# and the Cameron-Gelbach-Miller combination of three CR1 variances with two.

#' Cluster ids (1..G) of the rows lm() used, one vector per cluster
#' dimension, plus their intersection with two dimensions. NULL without
#' clusters.
#' @noRd
.se_groups <- function(fit, cluster_list, n_data) {
  if (is.null(cluster_list)) return(NULL)
  dropped <- as.integer(fit$na.action)
  keep    <- if (length(dropped)) setdiff(seq_len(n_data), dropped) else seq_len(n_data)
  g <- lapply(cluster_list, function(v) { v <- v[keep]; match(v, unique(v)) })
  if (length(g) == 2) {
    both   <- paste(g[[1]], g[[2]])
    g[[3]] <- match(both, unique(both))
  }
  g
}

#' Conventional standard error of the treatment coefficient and the pieces
#' its sandwich variances are built from: the score of each observation
#' (weight x residualised treatment x outcome residual) and the weighted sum
#' of squares of the residualised treatment.
#' @noRd
.se_core <- function(y, X, w, treat) {
  wfit <- function(X, y) if (is.null(w)) stats::lm.fit(X, y) else stats::lm.wfit(X, y, w)
  f    <- wfit(X, y)
  ww   <- if (is.null(w)) rep(1, length(y)) else w
  n    <- sum(ww > 0); k <- f$rank; df <- n - k
  dres <- wfit(X[, colnames(X) != treat, drop = FALSE], X[, treat])$residuals
  den  <- sum(ww * dres^2)
  list(se_conv = sqrt(sum(ww * f$residuals^2) / df / den),
       score = ww * dres * f$residuals, den = den, n = n, k = k, df = df)
}

#' Robust standard error of the treatment coefficient from .se_core():
#' HC1 (no clusters), CR1 (one cluster variable) or two-way CR1 by inclusion
#' and exclusion (two). `k_extra` adds regressors to the finite-sample
#' correction, as when an omitted variable is thought of as included. NA if
#' the two-way variance is not positive.
#' @noRd
.robust_se <- function(core, groups, k_extra = 0) {
  n <- core$n; p <- core$k + k_extra
  if (n - p <= 0) return(NA_real_)
  cr <- function(g) {
    G <- max(g)
    if (G < 2) return(NA_real_)
    sum(rowsum(core$score, g)^2) / core$den^2 * G / (G - 1) * (n - 1) / (n - p)
  }
  v <- if (is.null(groups)) sum(core$score^2) / core$den^2 * n / (n - p)
       else if (length(groups) == 1) cr(groups[[1]])
       else cr(groups[[1]]) + cr(groups[[2]]) - cr(groups[[3]])
  if (!is.finite(v) || v <= 0) NA_real_ else sqrt(v)
}

.robust_label <- function(groups, cluster_names) {
  if (is.null(groups)) return("HC1")
  base <- if (length(groups) == 1) "CR1" else "two-way CR1"
  if (is.null(cluster_names)) base else
    sprintf("%s by %s", base, paste(cluster_names, collapse = " and "))
}

#' sensemakr's adjusted estimate written in the coefficient, its scale and
#' the bias factor bf = sqrt(r2yz.dx * r2dz.x / (1 - r2dz.x)); `bf` may be a
#' vector, in which case a matrix with one column per value is returned.
#' @noRd
.adj_bf <- function(est, scale, bf, reduce = TRUE) {
  s <- if (isTRUE(reduce)) -1 else 1
  if (length(bf) == 1) sign(est) * (abs(est) + s * scale * bf)
  else sign(est) * (abs(est) + s * outer(scale, bf))
}

#' Bootstrap standard error of the adjusted estimate for an omitted variable
#' of fixed strength, from the stored coefficient and scale of every
#' replicate: what the table reports when the strength is not re-estimated.
#' @noRd
.fixed_strength_se <- function(boot_base, bf, reduce, is_twoway) {
  if (!is.finite(bf)) return(NA_real_)
  one <- function(m) {
    a <- .adj_bf(m[, 1], m[, 2], bf, reduce); a <- a[is.finite(a)]
    if (length(a) < 10) NA_real_ else var(a)
  }
  if (!is_twoway) return(sqrt(one(boot_base)))
  v <- vapply(boot_base, one, numeric(1))
  if (anyNA(v)) return(NA_real_)
  v2 <- v[["G"]] + v[["H"]] - v[["GH"]]
  if (!is.finite(v2) || v2 < 0) v2 <- max(v[["G"]], v[["H"]])
  sqrt(v2)
}

#' Standard errors of the adjusted estimate, one row per kd
#' @noRd
.se_compare <- function(fit, treat, results, groups, boot_base, reduce,
                        is_twoway, bench_reestimated, has_strata, cluster_names) {
  core <- .se_core(stats::model.response(stats::model.frame(fit)),
                   model.matrix(fit), fit$weights, treat)
  dof  <- fit$df.residual
  r2d  <- results$r2dz.x; r2y <- results$r2yz.dx
  ok   <- is.finite(r2d) & is.finite(r2y) & r2d < 1
  shrink <- ifelse(ok, sqrt((1 - r2y) / (1 - r2d)), NA_real_)
  bf     <- ifelse(ok, sqrt(r2y * r2d / (1 - r2d)), NA_real_)
  tab <- data.frame(
    kd              = results$kd,
    bootstrap       = results$se,
    bootstrap_fixed = vapply(bf, function(b)
      .fixed_strength_se(boot_base, b, reduce, is_twoway), numeric(1)),
    conventional    = shrink * core$se_conv * sqrt(dof / (dof - 1)),
    robust          = shrink * .robust_se(core, groups, k_extra = 1))
  notes <- c(
    if (has_strata) "The closed-form robust standard error ignores the strata.",
    if (length(groups) == 3) paste(
      "Wainstein and Hazlett (2026) give the closed form for one cluster",
      "variable; its use with two-way clustering is an extension made here."))
  list(table = tab,
       unadjusted = c(conventional = core$se_conv, robust = .robust_se(core, groups)),
       robust_type = .robust_label(groups, cluster_names),
       bench_reestimated = bench_reestimated,
       note = if (length(notes)) paste(notes, collapse = " ") else NULL)
}

#' Ratio of the robust to the conventional standard error of the treatment
#' coefficient, and how much it changes when each benchmark is taken out of
#' the regression and restored (Acedanski, 2026, Section 6.2). `cap` is the
#' ratio plus the largest positive change. Benchmarks are matched to columns
#' of the model matrix; a grouped benchmark is removed as one block.
#' @noRd
.se_ratio <- function(fit, treat, groups, benchmark_covariates,
                      gbenchmark_covariates, cluster_names) {
  y <- stats::model.response(stats::model.frame(fit))
  X <- model.matrix(fit); w <- fit$weights
  ratio <- function(X) {
    core <- .se_core(y, X, w, treat)
    c(.robust_se(core, groups), core$se_conv)
  }
  full <- ratio(X)
  out  <- list(omega = full[1] / full[2], se_conventional = full[2],
               se_robust = full[1],
               robust_type = .robust_label(groups, cluster_names),
               benchmarks = NULL, cap = NA_real_)
  pool <- if (!is.null(gbenchmark_covariates)) list(gbenchmark_covariates)
          else as.list(benchmark_covariates)
  if (length(pool) == 0) return(out)
  without <- vapply(pool, function(b) {
    keep <- !(colnames(X) %in% b)
    if (all(keep)) return(NA_real_)
    r <- ratio(X[, keep, drop = FALSE])
    r[1] / r[2]
  }, numeric(1))
  out$benchmarks <- data.frame(
    benchmark = vapply(pool, paste, character(1), collapse = " + "),
    omega_without = without, change = out$omega - without,
    stringsAsFactors = FALSE)
  ch <- out$benchmarks$change[is.finite(out$benchmarks$change)]
  if (length(ch) && is.finite(out$omega)) out$cap <- out$omega + max(0, ch)
  out
}

.print_se_compare <- function(sc, labels) {
  if (is.null(sc)) return(invisible(NULL))
  f <- function(v) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = 6))
  tab <- data.frame(Bootstrap = f(sc$table$bootstrap),
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (isTRUE(sc$bench_reestimated))
    tab[["Benchmark fixed"]] <- f(sc$table$bootstrap_fixed)
  tab[["Conventional"]] <- f(sc$table$conventional)
  tab[["Closed-form"]]  <- f(sc$table$robust)
  rownames(tab) <- make.unique(labels, sep = " #")
  cat(sprintf("\nStandard errors of the adjusted estimate (closed form: %s):\n",
              sc$robust_type))
  print(tab, right = TRUE, quote = FALSE)
  txt <- paste0(
    "(Bootstrap: as in the table above. ",
    if (isTRUE(sc$bench_reestimated)) paste0(
      "Benchmark fixed: the bootstrap with the strength of the omitted ",
      "variable held at its value in the original fit, computed from the ",
      "same replicates. "),
    "Conventional: sensemakr's adjusted standard error. Closed-form: the ",
    "robust standard error of the original fit, rescaled as in Wainstein and ",
    "Hazlett (2026) on the assumption that the omitted variable shrinks every ",
    "residual by the same proportion. The last two refer to the regression ",
    "that would include the omitted variable, the bootstrap to the sampling ",
    "variability of the adjusted estimate.)")
  cat(strwrap(txt, width = 78), sep = "\n")
  if (!is.null(sc$note)) cat(strwrap(paste("Note:", sc$note), width = 78), sep = "\n")
  invisible(NULL)
}

.print_se_ratio <- function(sr, treat) {
  if (is.null(sr) || !is.finite(sr$omega)) return(invisible(NULL))
  f <- function(v, flag = "") if (is.na(v)) "NA" else
    formatC(v, format = "f", digits = 4, flag = flag)
  cat(sprintf("\nRobust-to-conventional ratio of the standard error of %s (%s): %s\n",
              treat, sr$robust_type, f(sr$omega)))
  b <- sr$benchmarks
  if (!is.null(b)) for (i in seq_len(nrow(b)))
    cat(sprintf("  without %s: %s (change when it is restored: %s)\n",
                b$benchmark[i], f(b$omega_without[i]), f(b$change[i], "+")))
  cat(strwrap(paste0(
    "(Unadjusted regression. A ratio away from one means that robust and ",
    "conventional inference already differ before any omitted variable is ",
    "considered. The change shows how far restoring an observed covariate ",
    "moves the ratio, a yardstick for what an omitted variable might do; ",
    "see Acedanski, 2026.)"), width = 78), sep = "\n")
  invisible(NULL)
}


# ==============================================================================
# Fixed-strength curve and contour surface
# ==============================================================================
# For an omitted variable of given strength the adjusted estimate depends on
# the two partial R2s only through the bias factor
# bf = sqrt(r2yz.dx * r2dz.x / (1 - r2dz.x)). The bootstrap interval and
# p-value at that strength, held fixed across replicates, follow from the
# stored coefficient and scale of every replicate.

#' Adjusted estimate, bootstrap interval and p-value at each bias factor
#' @noRd
.strength_curve <- function(x, bf) {
  alpha  <- x$alpha
  reduce <- if (is.null(x$reduce)) TRUE else x$reduce
  df     <- if (is.null(x$crit_df)) Inf else x$crit_df
  est    <- as.numeric(.adj_bf(x$robustness$estimate, x$robustness$scale, bf, reduce))
  na     <- rep(NA_real_, length(bf))
  amat   <- function(m) {
    ok <- is.finite(m[, 1]) & is.finite(m[, 2])
    if (sum(ok) < 10) return(NULL)
    a <- .adj_bf(m[ok, 1], m[ok, 2], bf, reduce)
    if (is.matrix(a)) a else matrix(a, ncol = 1)
  }
  if (identical(x$method, "cgm_twoway")) {
    v <- lapply(x$boot_base, function(m) {
      a <- amat(m); if (is.null(a)) na else apply(a, 2, var) })
    v2  <- v$G + v$H - v$GH
    neg <- is.finite(v2) & v2 < 0
    v2[neg] <- pmax(v$G, v$H)[neg]
    se  <- sqrt(v2); cr <- .crit(alpha, df)
    return(data.frame(bf = bf, estimate = est, lower = est - cr * se,
                      upper = est + cr * se, pvalue = .pval2(est / se, df)))
  }
  a <- amat(x$boot_base)
  if (is.null(a))
    return(data.frame(bf = bf, estimate = est, lower = na, upper = na, pvalue = na))
  qs <- apply(a, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
  data.frame(bf = bf, estimate = est, lower = qs[1, ], upper = qs[2, ],
             pvalue = 2 * pmin(colMeans(a <= 0), colMeans(a >= 0)))
}

#' Surface over a grid of the two partial R2s: `value` is "estimate",
#' "pvalue" or "ci", the limit of the interval nearer zero.
#' @noRd
.contour_surface <- function(x, r2dz, r2yz, value = "ci") {
  bfg <- sqrt(outer(r2dz / (1 - r2dz), r2yz))
  bfs <- seq(0, max(bfg), length.out = 400)
  cv  <- .strength_curve(x, bfs)
  col <- switch(value, estimate = "estimate", pvalue = "pvalue",
                ci = if (isTRUE(cv$estimate[1] < 0)) "upper" else "lower")
  if (all(is.na(cv[[col]])))
    stop("Too few usable replicates for a contour plot.", call. = FALSE)
  z <- matrix(stats::approx(bfs, cv[[col]], xout = as.vector(bfg), rule = 2)$y,
              nrow = length(r2dz))
  list(r2dz.x = r2dz, r2yz.dx = r2yz, z = z, which = col)
}


# ==============================================================================
# Print method
# ==============================================================================

#' Print a bootmakr object
#'
#' Prints the table of bias-adjusted estimates with bootstrap standard
#' errors, confidence intervals and p-values (one row per \code{kd}),
#' followed by the descriptive benchmark-strength block, the breakdown kd,
#' the robustness values, the standard errors of the adjusted estimate
#' side by side, the ratio of the robust to the conventional standard
#' error and, if requested, the convergence diagnostics.
#'
#' @param x An object of class \code{"bootmakr"}, as returned by
#'   \code{\link{bootmakr}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{bootmakr}}, \code{\link{plot.bootmakr}}
#' @export
print.bootmakr <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)

  is_twoway <- identical(x$method, "cgm_twoway")

  cat(sprintf(
    "\nBootstrap sensitivity analysis (%s reps%s, n = %s",
    formatC(x$N_reps, big.mark = ",", format = "d"),
    if (is_twoway) " per dimension" else "",
    formatC(x$N, big.mark = ",", format = "d")
  ))
  if (!is.null(x$N_clust)) {
    if (is_twoway) {
      nm1 <- x$cluster_names[1] %||% "G"
      nm2 <- x$cluster_names[2] %||% "H"
      cat(sprintf(", two-way clusters: %s=%d, %s=%d, intersection=%d",
                  nm1, x$N_clust["G"], nm2, x$N_clust["H"], x$N_clust["GH"]))
    } else {
      cat(sprintf(", %d clusters", x$N_clust))
    }
  }
  cat(")\n")

  cat(sprintf("Benchmark: %s | kd = %s, ky = %s\n",
              x$bench_label, paste(x$kd, collapse = " "), paste(x$ky, collapse = " ")))
  if (isTRUE(x$fix_benchmark))
    cat("Benchmark strength held at its value in the original fit (fix_benchmark = TRUE)\n")
  if (isTRUE(x$weights_reestimated))
    cat("Weights re-estimated in every replication\n")
  if (isTRUE(x$semiweights))
    cat("Benchmark measured with semi-weights (Wainstein & Hazlett, 2026)\n")
  crit_df <- if (is.null(x$crit_df)) Inf else x$crit_df

  alpha  <- x$alpha
  ci_pct <- round((1 - alpha) * 100)
  res    <- x$results

  labels <- if (nrow(res) > 1) paste0(x$treat, " (kd=", res$kd, ")") else x$treat
  stars  <- ifelse(res$pvalue < 0.001, "***",
            ifelse(res$pvalue < 0.01,  "**",
            ifelse(res$pvalue < 0.05,  "*",
            ifelse(res$pvalue < 0.1,   ".", " "))))
  fmt_p  <- vapply(res$pvalue, function(p)
    if (is.na(p)) "NA" else if (p == 0) "0"
    else format.pval(p, digits = 3, eps = 2e-16),
    character(1))

  tab <- data.frame(
    Estimate   = formatC(res$estimate,  format = "f", digits = 6),
    `Std. Err` = formatC(res$se,        format = "f", digits = 6),
    lo         = formatC(res$ci_lower,  format = "f", digits = 6),
    hi         = formatC(res$ci_upper,  format = "f", digits = 6),
    `Pr(>|0|)` = fmt_p,
    ` `        = stars,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rownames(tab) <- labels
  colnames(tab)[3] <- paste0(alpha / 2 * 100, "%")
  colnames(tab)[4] <- paste0((1 - alpha / 2) * 100, "%")

  ci_label <- if (!is_twoway) "percentile" else if (is.finite(crit_df)) "t-approx"
              else "normal-approx"
  cat(sprintf("\nAdjusted estimates (%s %d%% CI):\n", ci_label, ci_pct))
  print(tab, right = TRUE, quote = FALSE)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  if (is_twoway) {
    if (is.finite(crit_df))
      cat(sprintf(paste("(H0: adjusted estimate = 0; SE = sqrt(V_G + V_H - V_GH),",
                        "CI = est +/- t*SE and p-value from\n a t distribution with",
                        "min(G, H) - 1 = %d degrees of freedom)\n"), as.integer(crit_df)))
    else
    cat("(H0: adjusted estimate = 0; SE = sqrt(V_G + V_H - V_GH),",
        "CI = est +/- z*SE, p = 2*(1 - Phi(|z|)))\n")
    if (any(res$var_neg_fix, na.rm = TRUE))
      cat("Note: V_G + V_H - V_GH was negative for one or more rows;",
          "replaced with max(V_G, V_H) (conservative).\n")
  } else {
    cat("(H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)\n")
  }

  if (x$N_fail > 0)
    cat(sprintf("\nNote: %d of %d replications failed and were dropped%s.\n",
                x$N_fail, if (is_twoway) 3L * x$N_reps else x$N_reps,
                if (is_twoway) " (summed across the 3 sub-bootstraps)" else ""))

  .print_benchmark_strength(x$benchmark_strength)

  .print_breakdown(x$breakdown)

  .print_robustness(x$robustness, crit_df)

  .print_se_compare(x$se_compare, labels)

  .print_se_ratio(x$se_ratio, x$treat)

  if (!is.null(x$convergence)) .print_convergence(x$convergence)
  invisible(x)
}

.print_convergence <- function(conv) {
  s <- conv$summary; opts <- conv$opts
  cat(sprintf("\nConvergence diagnostics (reps %d to %d by %d, threshold = %d):\n",
              opts$minreps, max(conv$data$reps), opts$stepsize, s$threshold))
  fmt <- function(x, d = "%.4f") if (is.na(x) || !is.finite(x)) "   ---" else sprintf(d, x)
  tab <- data.frame(
    Mean = c(sprintf("%.4f", s$se_mean),  sprintf("%.4f", s$p_mean)),
    Range = c(sprintf("%.4f", s$se_range), sprintf("%.4f", s$p_range)),
    `CV %` = c(sprintf("%.2f", s$se_cv),   sprintf("%.2f", s$p_cv)),
    `Range (hi)` = c(fmt(s$se_range_hi), fmt(s$p_range_hi)),
    `CV % (hi)`  = c(fmt(s$se_cv_hi, "%.2f"), fmt(s$p_cv_hi, "%.2f")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rownames(tab) <- c("Std. error", "P-value")
  colnames(tab)[4] <- sprintf("Range (>=%d)", s$threshold)
  colnames(tab)[5] <- sprintf("CV%% (>=%d)", s$threshold)
  print(tab, right = TRUE, quote = FALSE)
  cat("(CV % = coefficient of variation: sd / mean * 100)\n")
}


# ==============================================================================
# Plot method
# ==============================================================================

#' Plot a bootmakr object
#'
#' @param x An object of class \code{"bootmakr"}, as returned by
#'   \code{\link{bootmakr}}.
#' @param type Which plot to draw:
#'   \describe{
#'     \item{\code{"kd_sweep"}}{adjusted estimate and bootstrap confidence
#'       interval at each \code{kd}; solid markers are significant at
#'       \code{alpha}, hollow markers are not. Shows the "breakdown point"
#'       at which the estimate stops being statistically significant.}
#'     \item{\code{"convergence"}}{three panels: the bootstrap distribution,
#'       and the standard error and p-value as functions of the number of
#'       replications (requires \code{converge} in the call to
#'       \code{bootmakr()}).}
#'     \item{\code{"histogram"}}{the bootstrap distribution of the adjusted
#'       estimate with the original estimate marked.}
#'     \item{\code{"contour"}}{contours, over the two partial R2 of the
#'       omitted variable, of the limit of the bootstrap confidence interval
#'       that lies nearer zero, computed from the stored replicates with that
#'       strength held fixed across replicates. The dashed line marks where
#'       the bootstrap test stops rejecting; on the diagonal it passes through
#'       the bootstrap robustness value. The benchmark of the table is marked,
#'       but when the benchmark is re-estimated in every replicate the
#'       interval in the table differs from the one drawn here (see
#'       \code{fix_benchmark} in \code{\link{bootmakr}}). Not chosen by
#'       \code{"auto"}.}
#'     \item{\code{"auto"}}{(default) \code{"convergence"} if diagnostics
#'       are available, otherwise \code{"kd_sweep"} if several \code{kd}
#'       values were used, otherwise \code{"histogram"}.}
#'   }
#' @param ... Further arguments; for \code{type = "histogram"},
#'   \code{kd_idx} selects which \code{kd} value to show (default 1);
#'   for \code{type = "contour"}, \code{value} (\code{"ci"}, the default,
#'   \code{"estimate"} or \code{"pvalue"}), \code{lim} and \code{lim.y} (upper
#'   ends of the two axes), \code{grid} and \code{nlevels}.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{bootmakr}}, \code{\link{print.bootmakr}}
#' @export
plot.bootmakr <- function(x, type = c("auto", "kd_sweep", "convergence", "histogram",
                                      "contour"), ...) {
  type <- match.arg(type)
  if (type == "auto") {
    type <- if (!is.null(x$convergence)) "convergence"
            else if (length(x$kd) > 1) "kd_sweep"
            else "histogram"
  }
  switch(type,
    kd_sweep    = .plot_kd_sweep(x, ...),
    convergence = .plot_convergence(x, ...),
    histogram   = .plot_histogram(x, ...),
    contour     = .plot_contour(x, ...)
  )
  invisible(x)
}

.plot_contour <- function(x, value = c("ci", "estimate", "pvalue"),
                          lim = NULL, lim.y = NULL, grid = 100, nlevels = 10, ...) {
  value <- match.arg(value)
  if (is.null(x$boot_base) || is.null(x$robustness))
    stop("This object does not hold the replicates a contour plot needs.", call. = FALSE)
  res <- x$results
  pts <- is.finite(res$r2dz.x) & is.finite(res$r2yz.dx)
  if (is.null(lim))   lim   <- min(max(0.4, 1.2 * res$r2dz.x[pts]), 0.99)
  if (is.null(lim.y)) lim.y <- min(max(0.4, 1.2 * res$r2yz.dx[pts]), 1)
  if (!(lim > 0 && lim < 1 && lim.y > 0 && lim.y <= 1))
    stop("`lim` must lie between 0 and 1 (exclusive) and `lim.y` between 0 and 1.",
         call. = FALSE)
  s <- .contour_surface(x, seq(0, lim, length.out = grid),
                        seq(0, lim.y, length.out = grid), value)
  ci_pct <- round((1 - x$alpha) * 100)
  how    <- if (identical(x$method, "cgm_twoway")) "two-way bootstrap" else "percentile bootstrap"
  ttl <- switch(s$which,
    estimate = "Adjusted estimate",
    pvalue   = sprintf("P-value of the adjusted estimate (%s)", how),
    lower    = sprintf("Lower limit of the %d%% interval (%s)", ci_pct, how),
    upper    = sprintf("Upper limit of the %d%% interval (%s)", ci_pct, how))
  old_par <- par(mar = c(6.5, 5, 3.6, 1.5)); on.exit(par(old_par))
  xlab <- "Partial R2 of the omitted variable with the treatment"
  ylab <- "Partial R2 of the omitted variable with the outcome"
  zr   <- range(s$z, na.rm = TRUE)
  if (zr[1] == zr[2]) {
    # nothing to contour: the value is the same over the whole range
    plot(NA, xlim = c(0, lim), ylim = c(0, lim.y), xlab = xlab, ylab = ylab,
         main = ttl, las = 1, cex.lab = 1.05, cex.main = 1)
    graphics::text(lim / 2, lim.y / 2,
                   sprintf("Constant over this range: %s", format(zr[1], digits = 4)))
  } else {
    graphics::contour(s$r2dz.x, s$r2yz.dx, s$z, nlevels = nlevels, col = "gray35",
                      xlab = xlab, ylab = ylab, main = ttl, las = 1,
                      cex.lab = 1.05, cex.main = 1, ...)
  }
  # Where the test of the table stops rejecting: the bias factor of the
  # bootstrap robustness value, traced over all pairs that share it
  rv <- x$robustness$rv_qa_boot
  if (isTRUE(is.finite(rv) && rv > 0 && rv < 1)) {
    bf0 <- rv / sqrt(1 - rv)
    xs  <- seq(lim / 2000, lim, length.out = 800)
    graphics::lines(xs, bf0^2 * (1 - xs) / xs, col = "red", lwd = 2, lty = 2)
    mtext(sprintf("Dashed line: from here on the bootstrap test is not significant at %s%%",
                  format(x$alpha * 100)), side = 3, line = 0.25, cex = 0.75, col = "red")
  }
  points(0, 0, pch = 17, cex = 1.1)
  if (any(pts)) {
    # label with the benchmark the table refers to; a manual bound otherwise
    bench <- x$breakdown$benchmark
    lab <- if (isTRUE(!is.na(bench))) sprintf("%sx %s", format(res$kd[pts]), bench)
           else if (identical(x$benchmark_strength$type, "manual")) rep(x$bench_label, sum(pts))
           else rep("manual R2", sum(pts))
    points(res$r2dz.x[pts], res$r2yz.dx[pts], pch = 23, bg = "navy", col = "navy", cex = 1.1)
    graphics::text(res$r2dz.x[pts], res$r2yz.dx[pts], lab, pos = 4, cex = 0.75, col = "navy")
  }
  note <- "Strength of the omitted variable held fixed across replications."
  if (isTRUE(x$se_compare$bench_reestimated))
    note <- paste(note, "The table re-estimates the benchmark, so its interval differs.")
  mtext(note, side = 1, line = 4.6, cex = 0.72)
  invisible(s)
}

.plot_kd_sweep <- function(x, ...) {
  res <- x$results; alpha <- x$alpha
  old_par <- par(mar = c(5.5, 5, 2, 1.5)); on.exit(par(old_par))
  kr <- diff(range(res$kd)); xlim <- range(res$kd) + c(-0.05, 0.05) * max(kr, 1)
  ylim <- range(c(res$ci_lower, res$ci_upper), na.rm = TRUE)
  ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
  plot(res$kd, res$estimate, type = "n", xlim = xlim, ylim = ylim,
       xlab = sprintf("Benchmark strength (kd x %s)", x$bench_label),
       ylab = "Adjusted Treatment Effect",
       las = 1, cex.lab = 1.1)
  abline(h = 0, col = "gray60", lty = 2)
  segments(res$kd, res$ci_lower, res$kd, res$ci_upper, col = "navy", lwd = 2.5)
  sig <- res$pvalue < alpha
  points(res$kd[sig],  res$estimate[sig],  pch = 16, col = "navy", cex = 1.6)
  points(res$kd[!sig], res$estimate[!sig], pch = 1,  col = "navy", cex = 1.6)
  mtext(sprintf("Note: %d%% CI from %s bootstrap reps. Solid = p < %.2f",
                round((1 - alpha) * 100),
                formatC(x$N_reps, big.mark = ",", format = "d"), alpha),
        side = 1, line = 4, cex = 0.75)
}

.plot_convergence <- function(x, ...) {
  if (is.null(x$convergence)) { message("No convergence data."); return(invisible(x)) }
  conv <- x$convergence; obs <- conv$obs_estimate
  old_par <- par(mfrow = c(3, 1), mar = c(4.5, 4, 2.5, 1.5)); on.exit(par(old_par))
  hist(conv$boot_vals, breaks = 40, col = adjustcolor("navy", 0.3),
       border = "navy", main = "Bootstrap Distribution",
       xlab = "", ylab = "", yaxt = "n", las = 1, cex.lab = 1.1)
  abline(v = obs, lty = 2, lwd = 2)
  plot(conv$data$reps, conv$data$se, type = "b", pch = 16, col = "navy",
       xlab = "", ylab = "", main = "SE Convergence",
       las = 1, cex.lab = 1.1)
  plot(conv$data$reps, conv$data$pvalue, type = "b", pch = 16, col = "maroon",
       xlab = "Number of Bootstrap Replications", ylab = "",
       main = "P-value Convergence", las = 1, cex.lab = 1.1)
}

.plot_histogram <- function(x, kd_idx = 1, ...) {
  obs <- x$results$estimate[kd_idx]
  is_twoway <- identical(x$method, "cgm_twoway")

  if (is_twoway) {
    nm1 <- x$cluster_names[1] %||% "G"
    nm2 <- x$cluster_names[2] %||% "H"
    panels <- list(
      list(vals = x$boot_samples$G[,  kd_idx],
           title = sprintf("Cluster: %s", nm1)),
      list(vals = x$boot_samples$H[,  kd_idx],
           title = sprintf("Cluster: %s", nm2)),
      list(vals = x$boot_samples$GH[, kd_idx],
           title = sprintf("Intersection: %s x %s", nm1, nm2))
    )
    old_par <- par(mfrow = c(3, 1), mar = c(3.5, 2, 2.5, 1.5)); on.exit(par(old_par))
    for (p in panels) {
      vals <- p$vals[is.finite(p$vals)]
      hist(vals, breaks = 40, col = adjustcolor("navy", 0.3), border = "navy",
           main = p$title, xlab = "Adjusted Estimate",
           ylab = "", yaxt = "n", las = 1, cex.lab = 1.0)
      abline(v = obs, lty = 2, lwd = 2)
    }
    mtext(sprintf("Two-way CGM bootstrap (kd = %s); dashed = original estimate",
                  x$kd[kd_idx]),
          side = 1, line = 2.2, cex = 0.75)
    return(invisible(NULL))
  }

  vals <- x$boot_samples[, kd_idx]; vals <- vals[is.finite(vals)]
  old_par <- par(mar = c(5, 2, 3, 1.5)); on.exit(par(old_par))
  hist(vals, breaks = 40, col = adjustcolor("navy", 0.3), border = "navy",
       main = sprintf("Bootstrap Distribution (kd = %s)", x$kd[kd_idx]),
       xlab = "Adjusted Estimate", ylab = "", yaxt = "n",
       las = 1, cex.lab = 1.1)
  abline(v = obs, lty = 2, lwd = 2)
  legend("topright", "Original estimate", lty = 2, lwd = 2, bty = "n")
}
