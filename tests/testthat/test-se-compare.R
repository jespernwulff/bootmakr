# Standard errors of the adjusted estimate side by side, and the ratio of the
# robust to the conventional standard error of the treatment coefficient.

data("firms", package = "bootmakr")

core_of <- function(fit, treat = "x")
  bootmakr:::.se_core(model.response(model.frame(fit)), model.matrix(fit),
                      fit$weights, treat)

test_that("the robust standard errors are those of the sandwich package", {
  skip_if_not_installed("sandwich")
  firms$w <- 0.5 + (firms$year - min(firms$year)) / diff(range(firms$year))
  ref <- function(v) sqrt(v["x", "x"])
  for (weighted in c(FALSE, TRUE)) {
    fit  <- if (weighted) lm(y ~ x + c, firms, weights = w) else lm(y ~ x + c, firms)
    core <- core_of(fit)
    expect_equal(core$se_conv, summary(fit)$coefficients["x", "Std. Error"], tolerance = 1e-10)
    expect_equal(bootmakr:::.robust_se(core, NULL),
                 ref(sandwich::vcovHC(fit, type = "HC1")), tolerance = 1e-10)
    g1 <- bootmakr:::.se_groups(fit, list(firms$firm), nrow(firms))
    expect_equal(bootmakr:::.robust_se(core, g1),
                 ref(sandwich::vcovCL(fit, cluster = firms$firm, type = "HC1")), tolerance = 1e-10)
    g2 <- bootmakr:::.se_groups(fit, list(firms$firm, firms$year), nrow(firms))
    expect_equal(bootmakr:::.robust_se(core, g2),
                 ref(sandwich::vcovCL(fit, cluster = data.frame(firms$firm, firms$year),
                                      type = "HC1", multi0 = FALSE)), tolerance = 1e-10)
  }
})

test_that("clusters are matched to the rows the regression used", {
  skip_if_not_installed("sandwich")
  d <- firms; d$c[c(2, 40, 300)] <- NA
  fit <- lm(y ~ x + c, d)
  g   <- bootmakr:::.se_groups(fit, list(d$firm), nrow(d))
  expect_length(g[[1]], nobs(fit))
  expect_equal(bootmakr:::.robust_se(core_of(fit), g),
               sqrt(sandwich::vcovCL(fit, cluster = d$firm[-c(2, 40, 300)], type = "HC1")["x", "x"]),
               tolerance = 1e-10)
})

test_that("the block has the bootstrap, sensemakr's and the closed-form standard error", {
  out <- bootmakr(y ~ x + c, firms, treat = "x", benchmark_covariates = "c",
                  kd = c(1, 1.5), cluster = "firm", reps = 60, seed = 3, progress = FALSE)
  tab <- out$se_compare$table
  expect_identical(names(tab), c("kd", "bootstrap", "bootstrap_fixed", "conventional", "robust"))
  expect_equal(tab$bootstrap, out$results$se)
  expect_true(out$se_compare$bench_reestimated)
  expect_identical(out$se_compare$robust_type, "CR1 by firm")

  fit <- lm(y ~ x + c, firms)
  ob  <- sensemakr::ovb_bounds(fit, treatment = "x", benchmark_covariates = "c", kd = c(1, 1.5))
  expect_equal(tab$conventional, ob$adjusted_se, tolerance = 1e-10)

  # closed form: the CR1 standard error of the fit, rescaled, with one more
  # regressor in the finite-sample correction (Wainstein & Hazlett, 2026)
  n   <- nobs(fit)
  cr1 <- bootmakr:::.robust_se(core_of(fit), bootmakr:::.se_groups(fit, list(firms$firm), nrow(firms)))
  expect_equal(out$se_compare$unadjusted[["robust"]], cr1)
  expect_equal(tab$robust,
               sqrt((1 - ob$r2yz.dx) / (1 - ob$r2dz.x)) * cr1 * sqrt((n - 3) / (n - 4)),
               tolerance = 1e-10)

  expect_output(print(out), "Standard errors of the adjusted estimate")
  expect_output(print(out), "Benchmark fixed")
  expect_output(print(out), "(closed form: CR1 by firm)", fixed = TRUE)
})

test_that("labels and notes follow the resampling scheme", {
  hc <- bootmakr(y ~ x + c, firms, treat = "x", r2dz.x = 0.02, r2yz.dx = 0.05,
                 strata = "year", reps = 20, seed = 1, progress = FALSE)
  expect_identical(hc$se_compare$robust_type, "HC1")
  expect_match(hc$se_compare$note, "ignores the strata")
  expect_false(hc$se_compare$bench_reestimated)
  # with a manual R2 nothing is re-estimated: the two bootstrap columns agree
  expect_equal(hc$se_compare$table$bootstrap_fixed, hc$results$se, tolerance = 1e-10)
  expect_null(hc$se_ratio$benchmarks)

  tw <- suppressWarnings(bootmakr(y ~ x + c, firms, treat = "x", benchmark_covariates = "c",
                                  cluster = c("firm", "year"), reps = 20, seed = 1, progress = FALSE))
  expect_identical(tw$se_compare$robust_type, "two-way CR1 by firm and year")
  expect_match(tw$se_compare$note, "extension made here")
})

test_that("the ratio of robust to conventional standard error, and its change", {
  out <- bootmakr(y ~ x + c, firms, treat = "x", benchmark_covariates = "c",
                  cluster = "firm", reps = 20, seed = 1, progress = FALSE)
  sr  <- out$se_ratio
  expect_equal(sr$omega, sr$se_robust / sr$se_conventional)
  expect_equal(sr$se_conventional,
               summary(lm(y ~ x + c, firms))$coefficients["x", "Std. Error"], tolerance = 1e-10)

  fit0  <- lm(y ~ x, firms)
  core0 <- core_of(fit0)
  om0   <- bootmakr:::.robust_se(core0, bootmakr:::.se_groups(fit0, list(firms$firm), nrow(firms))) /
    core0$se_conv
  expect_identical(sr$benchmarks$benchmark, "c")
  expect_equal(sr$benchmarks$omega_without, om0, tolerance = 1e-10)
  expect_equal(sr$benchmarks$change, sr$omega - om0, tolerance = 1e-10)
  expect_equal(sr$cap, sr$omega + max(0, sr$benchmarks$change))
  expect_output(print(out), "Robust-to-conventional ratio")
  expect_output(print(out), "without c:")
})

test_that("a grouped benchmark is taken out as one block", {
  skip_helper()
  f   <- peacefactor ~ directlyharmed + age + pastvoted + female
  out <- bootmakr(f, darfur_data, treat = "directlyharmed",
                  gbenchmark_covariates = c("female", "pastvoted"),
                  reps = 20, seed = 1, progress = FALSE)
  sr  <- out$se_ratio
  expect_identical(sr$benchmarks$benchmark, "female + pastvoted")
  fit0  <- lm(peacefactor ~ directlyharmed + age, darfur_data)
  core0 <- core_of(fit0, "directlyharmed")
  expect_equal(sr$benchmarks$omega_without,
               bootmakr:::.robust_se(core0, NULL) / core0$se_conv, tolerance = 1e-10)
})

test_that("the strength columns of the results describe the row the table reports", {
  skip_helper()
  f    <- peacefactor ~ directlyharmed + age + pastvoted + female
  cols <- c("estimate", "r2dz.x", "r2yz.dx", "impact")
  one  <- bootmakr(f, darfur_data, treat = "directlyharmed", benchmark_covariates = "female",
                   reps = 20, seed = 1, progress = FALSE)
  # second benchmark selected with bounds_row
  two  <- bootmakr(f, darfur_data, treat = "directlyharmed",
                   benchmark_covariates = c("age", "female"), bounds_row = 2,
                   reps = 20, seed = 1, progress = FALSE)
  expect_equal(two$results[, cols], one$results[, cols], tolerance = 1e-12)
  # a manual bound next to a benchmark: the table reports the manual row
  man  <- bootmakr(f, darfur_data, treat = "directlyharmed", r2dz.x = 0.02, r2yz.dx = 0.05,
                   reps = 20, seed = 1, progress = FALSE)
  expect_warning(
    both <- bootmakr(f, darfur_data, treat = "directlyharmed", benchmark_covariates = "female",
                     r2dz.x = 0.02, r2yz.dx = 0.05, reps = 20, seed = 1, progress = FALSE),
    "Multiple benchmark sources")
  expect_equal(both$results[, cols], man$results[, cols], tolerance = 1e-12)
})
