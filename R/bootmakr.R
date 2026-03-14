#' Bootstrap Inference for Sensemakr Sensitivity Analysis
#'
#' Wraps \code{\link[sensemakr]{sensemakr}} in a bootstrap loop to produce
#' bootstrap standard errors, percentile confidence intervals, and p-values
#' for the bias-adjusted treatment effect.
#'
#' @param formula A formula for the OLS regression.
#' @param data A data frame.
#' @param treat Character: treatment variable name.
#' @param benchmark_covariates Character vector: individual benchmark covariates.
#' @param gbenchmark_covariates Character vector: grouped benchmark covariates.
#'   Uses \code{group_partial_r2} to compute joint R2 values.
#' @param kd,ky Numeric vectors of benchmark multipliers. \code{ky} defaults to \code{kd}.
#' @param q,alpha Numeric: proportion of effect / significance level.
#' @param r2dz.x,r2yz.dx Optional manual R2 values (skip benchmark computation).
#' @param bound_label Character label for bounds table.
#' @param reduce Logical: bias-reducing confounders (default TRUE).
#' @param bounds_row Integer: which row of bounds table (default 1).
#' @param reps Integer: bootstrap replications (default 1000).
#' @param seed Integer or NULL.
#' @param cluster,strata,weights Column name or vector.
#' @param dots Deprecated; use \code{progress} instead.
#' @param progress Logical: show a progress bar (default TRUE).
#' @param converge TRUE, FALSE, or list(minreps, stepsize, threshold).
#' @param verbose Logical.
#' @return Object of class \code{"bootmakr"}.
#'
#' @importFrom sensemakr sensemakr group_partial_r2 ovb_partial_r2_bound adjusted_estimate
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline hist legend mtext par plot points segments
#' @importFrom stats complete.cases lm nobs quantile reformulate sd terms
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
  cluster_vec <- .resolve_var(cluster, data, "cluster")
  strata_vec  <- .resolve_var(strata, data, "strata")
  weight_vec  <- .resolve_var(weights, data, "weights")

  conv_opts <- .parse_converge(converge, reps)
  if (!is.null(seed)) set.seed(seed)

  # ---- Fit original model ----
  fit_orig <- if (is.null(weight_vec)) {
    lm(formula, data = data)
  } else {
    lm(formula, data = data, weights = weight_vec)
  }
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

  # ---- Set up resampling ----
  resample_info <- .setup_resampling(data, cluster_vec, strata_vec)

  # ---- Bootstrap loop ----
  boot_mat <- matrix(NA_real_, nrow = reps, ncol = n_kd)
  colnames(boot_mat) <- paste0("kd_", kd)
  n_fail <- 0L

  # Progress bar
  show_progress <- isTRUE(progress)
  if (show_progress) {
    clust_msg <- if (!is.null(resample_info$n_clust))
      sprintf(", %d clusters", resample_info$n_clust) else ""
    cat(sprintf("Bootstrapping (%s reps%s)\n", formatC(reps, big.mark = ","), clust_msg))
    pb <- txtProgressBar(min = 0, max = reps, style = 3, width = 50)
  }

  for (b in seq_len(reps)) {
    if (show_progress) setTxtProgressBar(pb, b)

    boot_idx <- .resample_once(resample_info)
    d_boot   <- data[boot_idx, , drop = FALSE]

    boot_mat[b, ] <- tryCatch({
      fit_b <- if (is.null(weight_vec)) {
        lm(formula, data = d_boot)
      } else {
        lm(formula, data = d_boot, weights = weight_vec[boot_idx])
      }
      .get_adjusted_estimates(
        fit_b, d_boot, formula, treat, weight_vec[boot_idx],
        benchmark_covariates, gbenchmark_covariates,
        kd, ky, q, alpha, r2dz.x, r2yz.dx, bound_label, reduce, bounds_row
      )
    }, error = function(e) rep(NA_real_, n_kd))

    if (any(!is.finite(boot_mat[b, ]))) n_fail <- n_fail + 1L
  }
  if (show_progress) { close(pb); cat("\n") }

  # ---- Statistics ----
  results      <- .compute_boot_stats(boot_mat, obs_estimates, kd, alpha)
  n_successful <- sum(complete.cases(boot_mat))

  conv_out <- NULL
  if (conv_opts$do_converge) {
    conv_out <- .convergence_diagnostics(boot_mat[, 1], obs_estimates[1], conv_opts)
  }

  structure(
    list(
      results       = results,
      boot_samples  = boot_mat,
      convergence   = conv_out,
      call          = cl,
      N             = N,
      N_reps        = reps,
      N_successful  = n_successful,
      N_fail        = n_fail,
      N_clust       = resample_info$n_clust,
      kd            = kd,
      ky            = ky,
      alpha         = alpha,
      treat         = treat,
      benchmark_covariates  = benchmark_covariates,
      gbenchmark_covariates = gbenchmark_covariates,
      bench_label   = bench_label,
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
    # 1. Partial R2 of Y with Z_group given D, X  (from the outcome model)
    r2yxj_base <- sensemakr::group_partial_r2(fit, covariates = gbenchmark_covariates)

    # 2. Partial R2 of D with Z_group given X  (from a treatment model)
    #    Build formula: treat ~ all other RHS variables
    rhs_vars <- attr(terms(formula), "term.labels")
    rhs_no_treat <- setdiff(rhs_vars, treat)
    treat_formula <- reformulate(rhs_no_treat, response = treat)
    fit_d <- if (is.null(weight_vec)) {
      lm(treat_formula, data = data)
    } else {
      lm(treat_formula, data = data, weights = weight_vec)
    }
    r2dxj_base <- sensemakr::group_partial_r2(fit_d, covariates = gbenchmark_covariates)

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
# Helpers
# ==============================================================================

.resolve_var <- function(x, data, label) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1 && x %in% names(data)) return(data[[x]])
  if (length(x) == nrow(data)) return(x)
  stop(sprintf("`%s` must be a column name in data or vector of length nrow(data).", label))
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


# ==============================================================================
# Print method
# ==============================================================================

#' @export
print.bootmakr <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)

  cat(sprintf(
    "\nBootstrap sensitivity analysis (%s reps, n = %s",
    formatC(x$N_reps, big.mark = ","), formatC(x$N, big.mark = ",")
  ))
  if (!is.null(x$N_clust)) cat(sprintf(", %d clusters", x$N_clust))
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
    if (p == 0) "0" else format.pval(p, digits = 3, eps = 2e-16),
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

  cat(sprintf("\nAdjusted estimates (percentile %d%% CI):\n", ci_pct))
  print(tab, right = TRUE, quote = FALSE)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat("(H0: adjusted estimate = 0; CI and p-value from percentile bootstrap)\n")

  if (x$N_fail > 0)
    cat(sprintf("\nNote: %d of %d replications failed and were dropped.\n",
                x$N_fail, x$N_reps))

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
  old_par <- par(mfrow = c(3, 1), mar = c(4.5, 5, 2.5, 1.5)); on.exit(par(old_par))
  hist(conv$boot_vals, breaks = 40, col = adjustcolor("navy", 0.3),
       border = "navy", main = "Bootstrap Distribution",
       xlab = "", ylab = "Frequency", las = 1, cex.lab = 1.1)
  abline(v = obs, lty = 2, lwd = 2)
  plot(conv$data$reps, conv$data$se, type = "b", pch = 16, col = "navy",
       xlab = "", ylab = "Standard Error", main = "SE Convergence",
       las = 1, cex.lab = 1.1)
  plot(conv$data$reps, conv$data$pvalue, type = "b", pch = 16, col = "maroon",
       xlab = "Number of Bootstrap Replications", ylab = "Two-sided P-value",
       main = "P-value Convergence", las = 1, cex.lab = 1.1)
}

.plot_histogram <- function(x, kd_idx = 1, ...) {
  vals <- x$boot_samples[, kd_idx]; vals <- vals[is.finite(vals)]
  obs <- x$results$estimate[kd_idx]
  old_par <- par(mar = c(5, 5, 3, 1.5)); on.exit(par(old_par))
  hist(vals, breaks = 40, col = adjustcolor("navy", 0.3), border = "navy",
       main = sprintf("Bootstrap Distribution (kd = %s)", x$kd[kd_idx]),
       xlab = "Adjusted Estimate", ylab = "Frequency",
       las = 1, cex.lab = 1.1)
  abline(v = obs, lty = 2, lwd = 2)
  legend("topright", "Original estimate", lty = 2, lwd = 2, bty = "n")
}
