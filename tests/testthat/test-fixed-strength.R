# fix_benchmark, the fixed-strength curve behind the contour plot, t critical
# values for two-way clustering, and weights given as a function.

data("firms", package = "bootmakr")
firms$w <- 0.5 + (firms$year - min(firms$year)) / diff(range(firms$year))
cols <- c("estimate", "se", "ci_lower", "ci_upper", "pvalue")
run  <- function(...) suppressWarnings(bootmakr(y ~ x + c, firms, treat = "x", ..., progress = FALSE))

test_that("fix_benchmark equals a manual-R2 run at the bound of the original fit", {
  d0 <- run(benchmark_covariates = "c", cluster = "firm", reps = 60, seed = 5)
  fx <- run(benchmark_covariates = "c", cluster = "firm", reps = 60, seed = 5, fix_benchmark = TRUE)
  mn <- run(r2dz.x = d0$results$r2dz.x, r2yz.dx = d0$results$r2yz.dx,
            cluster = "firm", reps = 60, seed = 5)
  expect_equal(fx$results[, cols], mn$results[, cols], tolerance = 1e-10)
  # the companion column of the default run is the table of the fixed run
  expect_equal(d0$se_compare$table$bootstrap_fixed, fx$results$se, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(d0$results$se, fx$results$se)))
  expect_true(fx$fix_benchmark); expect_false(d0$fix_benchmark)
  expect_false(fx$se_compare$bench_reestimated)
  # same resamples, same coefficients; the benchmark draws are constant
  expect_equal(fx$boot_base, d0$boot_base)
  expect_true(all(fx$boot_bench[, 1] == fx$boot_bench[1, 1]))
  expect_output(print(fx), "Benchmark strength held at its value in the original fit")
  expect_error(run(benchmark_covariates = "c", fix_benchmark = NA, reps = 20), "TRUE or FALSE")
})

test_that("with fix_benchmark the breakdown kd refers to the fixed strength", {
  args <- list(benchmark_covariates = "c", cluster = "firm", reps = 100, seed = 5, fix_benchmark = TRUE)
  fx <- do.call(run, args)
  k  <- fx$breakdown$kd
  expect_true(is.finite(k) && k > 0)
  expect_equal(fx$breakdown$n_used, 100L)
  p <- function(m) do.call(run, c(args, list(kd = k * m)))$results$pvalue
  expect_lt(p(1 - 1e-7), 0.05)
  expect_gte(p(1 + 1e-7), 0.05)
})

test_that("the fixed-strength curve reproduces the fixed run and the robustness value", {
  bf <- function(r) sqrt(r$r2yz.dx * r$r2dz.x / (1 - r$r2dz.x))
  for (cl in list("firm", c("firm", "year"))) {
    d0 <- run(benchmark_covariates = "c", cluster = cl, reps = 60, seed = 9)
    fx <- run(benchmark_covariates = "c", cluster = cl, reps = 60, seed = 9, fix_benchmark = TRUE)
    cv <- bootmakr:::.strength_curve(d0, bf(d0$results))
    expect_equal(unlist(cv[, c("estimate", "lower", "upper", "pvalue")], use.names = FALSE),
                 unlist(fx$results[, cols[-2]], use.names = FALSE), tolerance = 1e-10)
    rv <- d0$robustness$rv_qa_boot
    r  <- rv * c(1 - 1e-6, 1 + 1e-6)
    pp <- bootmakr:::.strength_curve(d0, r / sqrt(1 - r))$pvalue
    expect_lt(pp[1], 0.05); expect_gte(pp[2], 0.05)
  }
})

test_that("the contour surface and plot", {
  d0 <- run(benchmark_covariates = "c", kd = c(0.5, 1), cluster = "firm", reps = 60, seed = 9)
  s  <- bootmakr:::.contour_surface(d0, seq(0, 0.4, length.out = 30), seq(0, 0.4, length.out = 25), "ci")
  expect_equal(dim(s$z), c(30L, 25L))
  expect_identical(s$which, "lower")              # positive estimate: lower limit
  expect_equal(s$z[1, 1], unname(quantile(d0$boot_base[, "estimate"], 0.025)), tolerance = 1e-10)
  expect_true(all(diff(s$z[, 25]) <= 1e-12))      # stronger omitted variable, lower limit
  direct <- bootmakr:::.strength_curve(d0, sqrt(s$r2yz.dx[20] * s$r2dz.x[17] / (1 - s$r2dz.x[17])))$lower
  expect_lt(abs(s$z[17, 20] - direct), 1e-4)

  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  for (v in c("ci", "estimate", "pvalue"))
    expect_error(plot(d0, type = "contour", value = v), NA)
  expect_error(plot(d0, type = "contour", lim = 1.2), "lim")
  nr <- run(benchmark_covariates = "c", reduce = FALSE, reps = 30, seed = 1)
  expect_error(plot(nr, type = "contour"), NA)
})

test_that("twoway_crit = 't' reaches the table, the robustness value and the breakdown kd", {
  args <- list(cluster = c("firm", "year"), reps = 60, seed = 2)
  z  <- do.call(run, c(args, list(benchmark_covariates = "c")))
  tt <- do.call(run, c(args, list(benchmark_covariates = "c", twoway_crit = "t")))
  df <- min(z$N_clust[c("G", "H")]) - 1
  expect_true(is.infinite(z$crit_df))
  expect_equal(tt$crit_df, df)
  r <- tt$results
  expect_equal(r$se, z$results$se)
  expect_equal(r$ci_lower, r$estimate - qt(0.975, df) * r$se)
  expect_equal(r$ci_upper, r$estimate + qt(0.975, df) * r$se)
  expect_equal(r$pvalue, 2 * pt(-abs(r$estimate / r$se), df))
  expect_gt(r$pvalue, z$results$pvalue)

  rv <- tt$robustness$rv_qa_boot
  expect_lt(rv, z$robustness$rv_qa_boot)
  p_rv <- function(m) do.call(run, c(args, list(r2dz.x = rv * m, r2yz.dx = rv * m,
                                                 twoway_crit = "t")))$results$pvalue
  expect_lt(p_rv(1 - 1e-6), 0.05); expect_gte(p_rv(1 + 1e-6), 0.05)

  k <- tt$breakdown$kd
  expect_lt(k, z$breakdown$kd)
  p_kd <- function(m) do.call(run, c(args, list(benchmark_covariates = "c", kd = k * m,
                                                 twoway_crit = "t")))$results$pvalue
  expect_lt(p_kd(1 - 1e-6), 0.05); expect_gte(p_kd(1 + 1e-6), 0.05)

  expect_output(print(tt), "t-approx")
  expect_output(print(tt), "two-way t-test")
  expect_output(print(z), "normal-approx")
  # no effect without two-way clustering
  one <- run(benchmark_covariates = "c", cluster = "firm", reps = 20, seed = 2, twoway_crit = "t")
  expect_true(is.infinite(one$crit_df))
  expect_error(run(benchmark_covariates = "c", twoway_crit = "chi", reps = 20))
})

test_that("a weight function that returns the stored weights equals fixed weights", {
  for (cl in list(NULL, "firm", c("firm", "year"))) {
    a <- run(benchmark_covariates = "c", cluster = cl, weights = "w", reps = 40, seed = 6)
    b <- run(benchmark_covariates = "c", cluster = cl, weights = function(d) d$w, reps = 40, seed = 6)
    expect_equal(b$results, a$results)
    expect_equal(b$boot_samples, a$boot_samples)
    expect_equal(b$boot_bench, a$boot_bench)
    expect_equal(b$breakdown, a$breakdown)
    expect_true(b$weights_reestimated); expect_false(a$weights_reestimated)
  }
  expect_output(print(b), "Weights re-estimated in every replication")
})

test_that("constant weights from a function give the unweighted results", {
  a <- run(benchmark_covariates = "c", cluster = "firm", reps = 40, seed = 6)
  b <- run(benchmark_covariates = "c", cluster = "firm", weights = function(d) rep(1, nrow(d)),
           reps = 40, seed = 6)
  expect_equal(b$results, a$results, tolerance = 1e-8)
  expect_equal(b$boot_samples, a$boot_samples, tolerance = 1e-8)
})

test_that("weights are re-estimated on the resample", {
  # weights that depend on the sample they are computed on, here through
  # the mean of x in the data the function is given (weights from year
  # shares would not do: the panel is balanced, so they are constant)
  wfun <- function(d) 1 / (1 + (d$x - mean(d$x))^2)
  seen <- integer(0)
  spy  <- function(d) { seen <<- c(seen, nrow(d)); wfun(d) }
  out  <- run(benchmark_covariates = "c", cluster = "firm", weights = spy, reps = 30, seed = 6)
  expect_length(seen, 31L)                      # original data and every resample
  fixed <- run(benchmark_covariates = "c", cluster = "firm", weights = wfun(firms), reps = 30, seed = 6)
  expect_equal(out$results$estimate, fixed$results$estimate)
  # the re-estimated weights reach the fit of every replicate
  expect_false(isTRUE(all.equal(out$boot_base, fixed$boot_base)))
  expect_false(isTRUE(all.equal(out$results$se, fixed$results$se)))
  expect_equal(out$N_fail, 0L)
})

test_that("a failing weight function is a failed replication; wrong output stops", {
  n <- 0
  flaky <- function(d) { n <<- n + 1; if (n %% 5 == 0) stop("no weights"); d$w }
  out <- run(benchmark_covariates = "c", weights = flaky, reps = 40, seed = 6)
  expect_equal(out$N_fail, 8L)
  expect_equal(sum(!is.finite(out$boot_base[, "estimate"])), 8L)
  msg <- "one non-negative number per row"
  expect_error(run(benchmark_covariates = "c", weights = function(d) 1, reps = 20), msg)
  expect_error(run(benchmark_covariates = "c", weights = function(d) rep(-1, nrow(d)), reps = 20), msg)
  expect_error(run(benchmark_covariates = "c", weights = function(d) letters[1], reps = 20), msg)
})
