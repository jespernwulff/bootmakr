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
#' omitted variable at each \code{kd}, together with their square roots
#' (absolute partial correlations). These are functions of the data alone.
#' No standard error enters them, so they do not depend on the variance
#' estimator. They are stored in \code{$benchmark_strength}.
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
#'   use a normal approximation \code{estimate +/- z * SE}. Two-way clustering
#'   is incompatible with \code{strata}.
#' @param strata,weights Column name or vector: strata for a stratified
#'   bootstrap, and regression weights.
#' @param dots Deprecated; use \code{progress} instead.
#' @param progress Logical: show a progress bar (default TRUE).
#' @param converge \code{TRUE}, \code{FALSE}, or
#'   \code{list(minreps, stepsize, threshold)}: an informal convergence check
#'   that recomputes the standard error and p-value on the first
#'   \code{minreps}, \code{minreps + stepsize}, ... replications and
#'   summarises how much they still move, overall and from \code{threshold}
#'   replications onwards.
#' @param verbose Logical; currently unused.
#' @return An object of class \code{"bootmakr"}, a list with (among others)
#'   \describe{
#'     \item{\code{results}}{data frame with one row per \code{kd}: adjusted
#'       \code{estimate}, bootstrap \code{se}, \code{ci_lower},
#'       \code{ci_upper}, \code{pvalue}, and the implied strength of the
#'       omitted variable (\code{r2dz.x}, \code{r2yz.dx}) the row is based
#'       on.}
#'     \item{\code{benchmark_strength}}{list with \code{$observed} (partial
#'       R2 of the benchmark with treatment and outcome, and their square
#'       roots) and \code{$implied} (sensemakr's bound on the omitted
#'       variable at each \code{kd}/\code{ky}, and the corresponding
#'       absolute partial correlations).}
#'     \item{\code{boot_samples}}{matrix of bootstrap draws of the adjusted
#'       estimate, one column per \code{kd} (a list of three such matrices
#'       with two-way clustering).}
#'     \item{\code{convergence}}{convergence diagnostics, if requested.}
#'     \item{\code{N}, \code{N_reps}, \code{N_successful}, \code{N_fail},
#'       \code{N_clust}}{sample size, replications requested, completed and
#'       failed, and number of clusters.}
#'     \item{\code{sensemakr_orig}}{the \code{sensemakr} object of the
#'       original fit.}
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
#' robust standard errors (and what to do instead). Working paper.
#'
#' @seealso \code{\link{print.bootmakr}}, \code{\link{plot.bootmakr}},
#'   \code{\link[sensemakr]{sensemakr}}
#'
#' @examples
#' data("darfur", package = "sensemakr")
#'
#' # An omitted variable as strong as `female` (kd = 1, the default);
#' # villages are resampled as whole clusters.
#' out <- bootmakr(
#'   peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
#'     pastvoted + hhsize_darfur + female + village,
#'   data = darfur, treat = "directlyharmed",
#'   benchmark_covariates = "female",
#'   cluster = "village", reps = 200, seed = 1, progress = FALSE
#' )
#' out
#'
#' # Benchmark strength on the partial-R2 and partial-correlation scales
#' out$benchmark_strength
#'
#' @importFrom sensemakr sensemakr group_partial_r2 ovb_partial_r2_bound adjusted_estimate partial_r2 ovb_bounds
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline hist legend mtext par plot points segments
#' @importFrom stats complete.cases lm model.matrix nobs pnorm qnorm quantile reformulate sd terms var
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
                     verbose = FALSE) {

  cl <- match.call()
  if (is.null(ky)) ky <- kd
  n_kd <- length(kd)

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
  weight_vec    <- .resolve_var(weights, data, "weights")
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
    kd, ky, q, alpha, r2dz.x, r2yz.dx, bound_label, reduce, bounds_row
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
    kd, ky, r2dz.x, r2yz.dx, bound_label
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

  # ---- Bootstrap loop(s) ----
  loop_args <- list(
    reps = reps, data = data, formula = formula, treat = treat,
    weight_vec = weight_vec,
    benchmark_covariates = benchmark_covariates,
    gbenchmark_covariates = gbenchmark_covariates,
    kd = kd, ky = ky, q = q, alpha = alpha,
    r2dz.x = r2dz.x, r2yz.dx = r2yz.dx,
    bound_label = bound_label, reduce = reduce, bounds_row = bounds_row,
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
    results      <- .compute_twoway_stats(boot_samples, obs_estimates, kd, alpha)
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
    results      <- .compute_boot_stats(boot_samples, obs_estimates, kd, alpha)
    n_fail       <- out$n_fail
    n_successful <- sum(complete.cases(boot_samples))
    n_clust      <- resample_info$n_clust
    method       <- "percentile"
  }
  results <- .attach_r2(results, bench_strength, n_kd)

  conv_out <- NULL
  if (conv_opts$do_converge) {
    conv_out <- if (is_twoway) {
      .convergence_diagnostics_twoway(boot_samples, obs_estimates[1], alpha, conv_opts)
    } else {
      .convergence_diagnostics(boot_samples[, 1], obs_estimates[1], conv_opts)
    }
  }

  structure(
    list(
      results       = results,
      boot_samples  = boot_samples,
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
                                    bound_label, reduce, bounds_row) {
  n_kd <- length(kd)

  if (!is.null(gbenchmark_covariates)) {
    # ---- Grouped benchmark pathway ----
    # 1.-2. Partial R2 of Y with Z_group given D, X (outcome model) and of
    #       D with Z_group given X (treatment model) -- see .group_base_r2()
    base <- .group_base_r2(fit, data, formula, treat, weight_vec, gbenchmark_covariates)
    r2yxj_base <- base$r2yxj.dx
    r2dxj_base <- base$r2dxj.x

    # 3. For each kd/ky, use ovb_partial_r2_bound for proper nonlinear scaling,
    #    then compute adjusted estimate with the scaled R2 values
    coefs  <- summary(fit)$coefficients
    est    <- coefs[treat, "Estimate"]
    se     <- coefs[treat, "Std. Error"]
    dof    <- fit$df.residual

    vals <- vapply(seq_len(n_kd), function(i) {
      bounds <- sensemakr::ovb_partial_r2_bound(
        r2dxj.x = r2dxj_base, r2yxj.dx = r2yxj_base,
        kd = kd[i], ky = ky[i], bound_label = "group"
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
  sm_args <- list(model = fit, treatment = treat)
  if (!is.null(benchmark_covariates)) sm_args$benchmark_covariates <- benchmark_covariates
  sm_args$kd    <- kd
  sm_args$ky    <- ky
  sm_args$q     <- q
  sm_args$alpha <- alpha
  if (!is.null(r2dz.x))     sm_args$r2dz.x     <- r2dz.x
  if (!is.null(r2yz.dx))    sm_args$r2yz.dx     <- r2yz.dx
  if (!is.null(bound_label)) sm_args$bound_label <- bound_label
  sm_args$reduce <- reduce

  sm <- do.call(sensemakr::sensemakr, sm_args)
  bnds <- sm$bounds
  if (is.null(bnds) || nrow(bnds) == 0) return(rep(NA_real_, n_kd))

  if (n_kd > 1) {
    bnds$adjusted_estimate[seq_len(n_kd)]
  } else {
    bnds$adjusted_estimate[bounds_row]
  }
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
#' descriptive block so that both report the same numbers.
#' @noRd
.group_base_r2 <- function(fit, data, formula, treat, weight_vec, gbenchmark_covariates) {
  r2yxj.dx <- sensemakr::group_partial_r2(fit, covariates = gbenchmark_covariates)
  rhs_vars     <- attr(terms(formula), "term.labels")
  rhs_no_treat <- setdiff(rhs_vars, treat)
  treat_formula <- reformulate(rhs_no_treat, response = treat)
  fit_d <- .fit_lm(treat_formula, data, weight_vec)
  r2dxj.x <- sensemakr::group_partial_r2(fit_d, covariates = gbenchmark_covariates)
  list(r2dxj.x = unname(r2dxj.x), r2yxj.dx = unname(r2yxj.dx))
}

#' Treatment model exactly as sensemakr's benchmarking routine builds it
#' (model matrix of the outcome model, treatment regressed on all other
#' columns, no intercept beyond the model-matrix constant, unweighted).
#' @noRd
.treatment_model_sm <- function(fit, treat) {
  m  <- model.matrix(fit)
  d  <- m[, treat]
  XX <- m[, !(colnames(m) %in% treat), drop = FALSE]
  lm(d ~ XX + 0)
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
#'              (sensemakr's bound on the omitted variable's strength) and
#'              their square roots r_dz.x, r_yz.dx.
#' All quantities are functions of the data alone (no standard errors), so
#' they are unaffected by heteroskedasticity- or cluster-robust inference.
#' Nothing here is a t-implied correlation or a threshold.
#' @noRd
.benchmark_strength <- function(fit, data, formula, treat, weight_vec,
                                benchmark_covariates, gbenchmark_covariates,
                                kd, ky, r2dz.x, r2yz.dx, bound_label) {
  n_kd <- length(kd)
  out  <- tryCatch({
    if (!is.null(gbenchmark_covariates)) {
      base  <- .group_base_r2(fit, data, formula, treat, weight_vec, gbenchmark_covariates)
      label <- paste(gbenchmark_covariates, collapse = " + ")
      observed <- data.frame(benchmark = label,
                             r2dxj.x = base$r2dxj.x, r2yxj.dx = base$r2yxj.dx,
                             stringsAsFactors = FALSE)
      implied <- do.call(rbind, lapply(seq_len(n_kd), function(i) {
        b <- sensemakr::ovb_partial_r2_bound(
          r2dxj.x = base$r2dxj.x, r2yxj.dx = base$r2yxj.dx,
          kd = kd[i], ky = ky[i], bound_label = "group"
        )
        data.frame(benchmark = label, kd = kd[i], ky = ky[i],
                   r2dz.x = b$r2dz.x, r2yz.dx = b$r2yz.dx,
                   stringsAsFactors = FALSE)
      }))
      list(type = "group", observed = observed, implied = implied)
    } else if (!is.null(benchmark_covariates)) {
      # Observed strength, computed as in sensemakr's benchmarking routine
      r2yxj.dx <- sensemakr::partial_r2(fit, covariates = benchmark_covariates)
      tm       <- .treatment_model_sm(fit, treat)
      r2dxj.x  <- sensemakr::partial_r2(tm, covariates = paste0("XX", benchmark_covariates))
      observed <- data.frame(benchmark = benchmark_covariates,
                             r2dxj.x = unname(r2dxj.x), r2yxj.dx = unname(r2yxj.dx),
                             stringsAsFactors = FALSE)
      # Implied strength of the omitted variable: sensemakr's own bounds
      # (rows ordered benchmark-by-benchmark, kd within benchmark -- the same
      # order used for the adjusted estimates)
      bnds <- sensemakr::ovb_bounds(fit, treatment = treat,
                                    benchmark_covariates = benchmark_covariates,
                                    kd = kd, ky = ky, adjusted_estimates = FALSE)
      n_b  <- length(benchmark_covariates)
      implied <- data.frame(benchmark = rep(benchmark_covariates, each = n_kd),
                            kd = rep(kd, times = n_b), ky = rep(ky, times = n_b),
                            r2dz.x = bnds$r2dz.x, r2yz.dx = bnds$r2yz.dx,
                            stringsAsFactors = FALSE)
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

  if (!is.null(out$implied)) {
    out$implied$r_dz.x  <- sqrt(out$implied$r2dz.x)
    out$implied$r_yz.dx <- sqrt(out$implied$r2yz.dx)
    rownames(out$implied) <- NULL
  }
  if (!is.null(out$observed)) {
    out$observed$r_dxj.x  <- sqrt(out$observed$r2dxj.x)
    out$observed$r_yxj.dx <- sqrt(out$observed$r2yxj.dx)
    rownames(out$observed) <- NULL
  }
  out
}

#' Attach the implied r2dz.x / r2yz.dx of each kd row to the results table
#' (first benchmark's rows, matching the adjusted estimates).
#' @noRd
.attach_r2 <- function(results, bs, n_kd) {
  results$r2dz.x  <- NA_real_
  results$r2yz.dx <- NA_real_
  imp <- bs$implied
  if (is.null(imp) || nrow(imp) == 0) return(results)
  if (identical(bs$type, "manual")) {
    results$r2dz.x  <- imp$r2dz.x[1]
    results$r2yz.dx <- imp$r2yz.dx[1]
  } else if (nrow(imp) >= n_kd) {
    results$r2dz.x  <- imp$r2dz.x[seq_len(n_kd)]
    results$r2yz.dx <- imp$r2yz.dx[seq_len(n_kd)]
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
  if (!is.null(ob) && nrow(ob) > 0) {
    rows <- unlist(lapply(seq_len(nrow(ob)), function(i) c(
      sprintf("%s with treatment | X", ob$benchmark[i]),
      sprintf("%s with outcome | D, X", ob$benchmark[i]))))
    r2   <- unlist(lapply(seq_len(nrow(ob)), function(i) c(ob$r2dxj.x[i], ob$r2yxj.dx[i])))
    tab  <- data.frame(`Partial R2` = f4(r2), `|Partial r|` = f3(sqrt(r2)),
                       check.names = FALSE, stringsAsFactors = FALSE)
    rownames(tab) <- rows
    print(tab, right = TRUE, quote = FALSE)
  }

  imp <- bs$implied
  if (!is.null(imp) && nrow(imp) > 0) {
    cat("Implied strength of the omitted variable (sensemakr bounds):\n")
    tab2 <- data.frame(kd = f2(imp$kd), ky = f2(imp$ky),
                       `R2dz.x` = f4(imp$r2dz.x), `R2yz.dx` = f4(imp$r2yz.dx),
                       `|r_dz.x|` = f3(imp$r_dz.x), `|r_yz.dx|` = f3(imp$r_yz.dx),
                       check.names = FALSE, stringsAsFactors = FALSE)
    rownames(tab2) <- if (identical(bs$type, "manual")) imp$benchmark else
      make.unique(sprintf("%sx %s", f2(imp$kd), imp$benchmark), sep = " #")
    print(tab2, right = TRUE, quote = FALSE)
  }
  cat("(R2dz.x, R2yz.dx: partial R2 of the omitted variable with the treatment given the\n",
      "covariates, and with the outcome given treatment and covariates; |r| = square root,\n",
      "i.e. the partial-correlation scale of the ITCV. Descriptive only: independent of the\n",
      "variance estimator. No t-implied correlation or threshold is reported.)\n", sep = " ")
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
                                bounds_row, show_progress, label = NULL) {
  n_kd     <- length(kd)
  boot_mat <- matrix(NA_real_, nrow = reps, ncol = n_kd)
  colnames(boot_mat) <- paste0("kd_", kd)
  n_fail   <- 0L

  if (show_progress) {
    clust_msg <- if (!is.null(resample_info$n_clust))
      sprintf(", %d clusters", resample_info$n_clust) else ""
    prefix <- if (!is.null(label)) paste0(label, " ") else ""
    cat(sprintf("%sBootstrapping (%s reps%s)\n", prefix,
                formatC(reps, big.mark = ","), clust_msg))
    pb <- txtProgressBar(min = 0, max = reps, style = 3, width = 50)
  }

  for (b in seq_len(reps)) {
    if (show_progress) setTxtProgressBar(pb, b)
    boot_idx <- .resample_once(resample_info)
    d_boot   <- data[boot_idx, , drop = FALSE]
    boot_mat[b, ] <- tryCatch({
      fit_b <- .fit_lm(formula, d_boot, weight_vec[boot_idx])
      .get_adjusted_estimates(
        fit_b, d_boot, formula, treat, weight_vec[boot_idx],
        benchmark_covariates, gbenchmark_covariates,
        kd, ky, q, alpha, r2dz.x, r2yz.dx, bound_label, reduce, bounds_row
      )
    }, error = function(e) rep(NA_real_, n_kd))
    if (any(!is.finite(boot_mat[b, ]))) n_fail <- n_fail + 1L
  }
  if (show_progress) { close(pb); cat("\n") }

  list(boot_mat = boot_mat, n_fail = n_fail)
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
.compute_twoway_stats <- function(boot_samples, obs_estimates, kd, alpha) {
  n_kd   <- length(kd)
  z_crit <- qnorm(1 - alpha / 2)
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
    results$pvalue[i]   <- 2 * pnorm(-abs(z_stat))
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
.convergence_diagnostics_twoway <- function(boot_samples, obs_estimate, alpha, opts) {
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
    conv_df$pvalue[j] <- 2 * pnorm(-abs(obs_estimate / se))
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
# Print method
# ==============================================================================

#' Print a bootmakr object
#'
#' Prints the table of bias-adjusted estimates with bootstrap standard
#' errors, confidence intervals and p-values (one row per \code{kd}),
#' followed by the descriptive benchmark-strength block and, if requested,
#' the convergence diagnostics.
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
    formatC(x$N_reps, big.mark = ","),
    if (is_twoway) " per dimension" else "",
    formatC(x$N, big.mark = ",")
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

  ci_label <- if (is_twoway) "normal-approx" else "percentile"
  cat(sprintf("\nAdjusted estimates (%s %d%% CI):\n", ci_label, ci_pct))
  print(tab, right = TRUE, quote = FALSE)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  if (is_twoway) {
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
#'     \item{\code{"auto"}}{(default) \code{"convergence"} if diagnostics
#'       are available, otherwise \code{"kd_sweep"} if several \code{kd}
#'       values were used, otherwise \code{"histogram"}.}
#'   }
#' @param ... Further arguments; for \code{type = "histogram"},
#'   \code{kd_idx} selects which \code{kd} value to show (default 1).
#' @return \code{x}, invisibly.
#' @seealso \code{\link{bootmakr}}, \code{\link{print.bootmakr}}
#' @export
plot.bootmakr <- function(x, type = c("auto", "kd_sweep", "convergence", "histogram"), ...) {
  type <- match.arg(type)
  if (type == "auto") {
    type <- if (!is.null(x$convergence)) "convergence"
            else if (length(x$kd) > 1) "kd_sweep"
            else "histogram"
  }
  switch(type,
    kd_sweep    = .plot_kd_sweep(x, ...),
    convergence = .plot_convergence(x, ...),
    histogram   = .plot_histogram(x, ...)
  )
  invisible(x)
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
                formatC(x$N_reps, big.mark = ","), alpha),
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
