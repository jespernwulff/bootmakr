# The bootstrap robustness value: the smallest strength r2dz.x = r2yz.dx at
# which the bootstrap test of the adjusted estimate stops rejecting. The
# defining property is checked by running bootmakr() again at that strength
# with the same seed.

data("firms", package = "bootmakr")

rerun_p <- function(r, ...) {
  bootmakr(..., r2dz.x = r, r2yz.dx = r, progress = FALSE)$results$pvalue
}

test_that("the robustness block has sensemakr's values and one row of base quantities per replicate", {
  skip_helper()
  out <- bootmakr(darfur_formula, data = darfur_data, treat = "directlyharmed",
                  benchmark_covariates = "female", cluster = "village",
                  reps = 80, seed = 3, progress = FALSE)
  rb  <- out$robustness
  fit <- lm(darfur_formula, data = darfur_data)

  expect_equal(rb$rv_q, as.numeric(sensemakr::robustness_value(
    fit, covariates = "directlyharmed")), tolerance = 1e-12)
  expect_equal(rb$rv_qa, as.numeric(sensemakr::robustness_value(
    fit, covariates = "directlyharmed", alpha = 0.05)), tolerance = 1e-12)
  expect_true(rb$rv_qa_boot > 0 && rb$rv_qa_boot < rb$rv_q)
  expect_null(rb$note)
  expect_equal(rb$n_used, 80L)

  expect_equal(dim(out$boot_base), c(80L, 2L))
  expect_identical(colnames(out$boot_base), c("estimate", "scale"))
  # the adjusted estimate of every replicate follows from the two columns
  bs <- out$benchmark_strength$implied
  b  <- out$boot_base
  manual <- bootmakr(darfur_formula, data = darfur_data, treat = "directlyharmed",
                     r2dz.x = 0.03, r2yz.dx = 0.03, cluster = "village",
                     reps = 80, seed = 3, progress = FALSE)
  expect_equal(unname(manual$boot_samples[, 1]),
               unname(sign(b[, 1]) * (abs(b[, 1]) - b[, 2] * 0.03 / sqrt(1 - 0.03))),
               tolerance = 1e-10)
  expect_true(is.data.frame(bs))
})

test_that("at the bootstrap robustness value the bootstrap p-value reaches alpha", {
  skip_helper()
  args <- list(darfur_formula, data = darfur_data, treat = "directlyharmed",
               cluster = "village", reps = 100, seed = 4)
  out <- do.call(bootmakr, c(args, list(benchmark_covariates = "female", progress = FALSE)))
  rv  <- out$robustness$rv_qa_boot

  expect_lt(do.call(rerun_p, c(list(rv * (1 - 1e-8)), args)), 0.05)
  expect_gte(do.call(rerun_p, c(list(rv * (1 + 1e-8)), args)), 0.05)
})

test_that("the same holds with weights and with another alpha", {
  skip_helper()
  d <- darfur_data
  d$w <- seq(0.5, 1.5, length.out = nrow(d))
  args <- list(peacefactor ~ directlyharmed + age + female + pastvoted, data = d,
               treat = "directlyharmed", weights = "w", alpha = 0.10,
               reps = 100, seed = 5)
  out <- do.call(bootmakr, c(args, list(gbenchmark_covariates = "female", progress = FALSE)))
  rv  <- out$robustness$rv_qa_boot

  expect_true(rv > 0)
  expect_lt(do.call(rerun_p, c(list(rv * (1 - 1e-8)), args)), 0.10)
  expect_gte(do.call(rerun_p, c(list(rv * (1 + 1e-8)), args)), 0.10)
})

test_that("with two-way clustering it is where the z-test stops rejecting", {
  args <- list(y ~ x + c, data = firms, treat = "x", cluster = c("firm", "year"),
               reps = 40, seed = 6)
  out <- do.call(bootmakr, c(args, list(benchmark_covariates = "c", progress = FALSE)))
  rv  <- out$robustness$rv_qa_boot

  expect_identical(out$robustness$method, "cgm_twoway")
  expect_true(is.list(out$boot_base) && identical(names(out$boot_base), c("G", "H", "GH")))
  expect_true(rv > 0 && rv < out$robustness$rv_q)
  expect_lt(do.call(rerun_p, c(list(rv * (1 - 1e-4)), args)), 0.05)
  expect_gt(do.call(rerun_p, c(list(rv * (1 + 1e-4)), args)), 0.05)
})

test_that("replicates with a coefficient of the opposite sign are counted on both sides", {
  # a weak effect: several replicates have a negative coefficient, so the
  # p-value is not monotone in the strength of the omitted variable
  set.seed(11)
  d <- data.frame(x = rnorm(150), c = rnorm(150))
  d$y <- 0.10 * d$x + 0.5 * d$c + rnorm(150)
  args <- list(y ~ x + c, data = d, treat = "x", alpha = 0.30, reps = 60, seed = 11)
  out <- do.call(bootmakr, c(args, list(benchmark_covariates = "c", progress = FALSE)))
  rv  <- out$robustness$rv_qa_boot

  expect_gt(sum(out$boot_base[, "estimate"] < 0), 0)
  expect_gt(rv, 0)
  expect_lt(do.call(rerun_p, c(list(rv * (1 - 1e-8)), args)), 0.30)
  expect_gte(do.call(rerun_p, c(list(rv * (1 + 1e-8)), args)), 0.30)

  # first strength on a grid at which the p-value reaches alpha
  b    <- out$boot_base
  grid <- seq(0, 0.2, length.out = 20001)
  p    <- vapply(grid, function(r) {
    a <- sign(b[, 1]) * (abs(b[, 1]) - b[, 2] * r / sqrt(1 - r))
    2 * min(mean(a <= 0), mean(a >= 0))
  }, numeric(1))
  expect_lt(abs(grid[which(p >= 0.30)[1]] - rv), 1.5e-5)
})

test_that("replicates whose benchmark bound fails still count", {
  # just below the kd at which the bound on r2dz.x reaches one, the bound is
  # impossible in part of the replicates; the robustness value does not use
  # the benchmark, so it is based on all of them and equals the value of a
  # run in which nothing fails
  r2c    <- unname(sensemakr::partial_r2(lm(x ~ c, data = firms), covariates = "c"))
  kd_max <- (1 - r2c) / r2c
  part <- suppressWarnings(bootmakr(y ~ x + c, data = firms, treat = "x",
                                    benchmark_covariates = "c", kd = 0.985 * kd_max,
                                    cluster = "firm", reps = 60, seed = 11, progress = FALSE))
  full <- bootmakr(y ~ x + c, data = firms, treat = "x", benchmark_covariates = "c",
                   cluster = "firm", reps = 60, seed = 11, progress = FALSE)

  expect_gt(part$N_fail, 0)
  expect_equal(part$robustness$n_used, 60L)
  expect_equal(part$robustness$rv_qa_boot, full$robustness$rv_qa_boot, tolerance = 1e-12)
})

test_that("cases without a bootstrap robustness value say why", {
  skip_helper()
  f <- peacefactor ~ directlyharmed + age + female
  half <- bootmakr(f, data = darfur_data, treat = "directlyharmed",
                   benchmark_covariates = "female", q = 0.5,
                   reps = 30, seed = 7, progress = FALSE)
  expect_true(is.na(half$robustness$rv_qa_boot))
  expect_match(half$robustness$note, "q = 1")

  up <- bootmakr(f, data = darfur_data, treat = "directlyharmed",
                 benchmark_covariates = "female", reduce = FALSE,
                 reps = 30, seed = 7, progress = FALSE)
  expect_true(is.na(up$robustness$rv_qa_boot))
  expect_match(up$robustness$note, "reduce = TRUE")

  # an estimate that is not significant to begin with has a value of zero
  set.seed(5)
  d <- data.frame(x = rnorm(150), c = rnorm(150))
  d$y <- 0.5 * d$c + rnorm(150)
  ns <- bootmakr(y ~ x + c, data = d, treat = "x", benchmark_covariates = "c",
                 reps = 60, seed = 8, progress = FALSE)
  expect_gte(rerun_p(0, y ~ x + c, data = d, treat = "x", reps = 60, seed = 8), 0.05)
  expect_identical(ns$robustness$rv_qa_boot, 0)
})

test_that("print shows the robustness values", {
  skip_helper()
  out <- bootmakr(peacefactor ~ directlyharmed + age + female, data = darfur_data,
                  treat = "directlyharmed", benchmark_covariates = "female",
                  reps = 30, seed = 9, progress = FALSE)
  txt <- capture.output(print(out))
  expect_true(any(grepl("Robustness value", txt, fixed = TRUE)))
  expect_true(any(grepl("Not significant at 5%, bootstrap", txt, fixed = TRUE)))
})
