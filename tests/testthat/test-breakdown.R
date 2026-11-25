# The breakdown kd: the smallest kd (with ky = kd) at which the bootstrap test
# of the adjusted estimate stops rejecting. The defining property is checked
# by running bootmakr() again at that kd with the same seed.

data("firms", package = "bootmakr")

p_at <- function(kd, ...) {
  suppressWarnings(bootmakr(..., kd = kd, progress = FALSE))$results$pvalue
}

# a benchmark strongly related to the treatment: sensemakr's bound exists only
# up to a kd of about 1.1, and less than that in many replicates
coll_data <- function() {
  set.seed(2)
  c <- rnorm(60); x <- 0.9 * c + rnorm(60); y <- 8 * x + 0.05 * c + rnorm(60)
  data.frame(y, x, c)
}

test_that("the closed-form bound is sensemakr's", {
  for (r2d in c(0.003, 0.14, 0.45)) for (r2y in c(0.001, 0.12, 0.5)) for (k in c(0.3, 1, 2.5)) {
    sm <- suppressWarnings(tryCatch(
      sensemakr::ovb_partial_r2_bound(r2dxj.x = r2d, r2yxj.dx = r2y, kd = k),
      error = function(e) NULL))
    me <- bootmakr:::.bound_at(r2d, r2y, k)
    if (is.null(sm)) {
      expect_true(is.na(me$r2dz.x) && is.na(me$r2yz.dx))
    } else {
      expect_equal(c(me$r2dz.x, me$r2yz.dx), c(sm$r2dz.x, sm$r2yz.dx), tolerance = 1e-12)
    }
  }
})

test_that("the breakdown block is complete and the bootstrap p-value reaches alpha there", {
  skip_helper()
  f    <- peacefactor ~ directlyharmed + age + pastvoted + female
  args <- list(f, data = darfur_data, treat = "directlyharmed",
               benchmark_covariates = "female", reps = 100, seed = 4)
  out <- do.call(bootmakr, c(args, list(progress = FALSE)))
  bd  <- out$breakdown

  expect_identical(bd$benchmark, "female")
  expect_identical(bd$method, "percentile")
  expect_equal(bd$n_used, 100L)
  expect_null(bd$note)
  expect_true(bd$kd > 0 && bd$kd < bd$kd_limit)
  expect_equal(dim(out$boot_bench), c(100L, 2L))
  expect_identical(colnames(out$boot_bench), c("r2dxj.x", "r2yxj.dx"))
  expect_true(all(out$boot_bench > 0 & out$boot_bench < 1))
  # $boot_base keeps its two columns
  expect_identical(colnames(out$boot_base), c("estimate", "scale"))

  expect_lt(do.call(p_at, c(list(bd$kd * (1 - 1e-8)), args)), 0.05)
  expect_gte(do.call(p_at, c(list(bd$kd * (1 + 1e-8)), args)), 0.05)

  # every kd below it is significant, in one run over several kd
  sweep <- do.call(bootmakr, c(args, list(kd = bd$kd * c(0.1, 0.5, 0.9, 0.99), progress = FALSE)))
  expect_true(all(sweep$results$pvalue < 0.05))
  # a kd sweep reports the same breakdown as a single kd
  expect_equal(sweep$breakdown$kd, bd$kd, tolerance = 1e-12)
})

test_that("the conventional breakdown kd is where sensemakr's interval reaches zero", {
  skip_helper()
  f   <- peacefactor ~ directlyharmed + age + pastvoted + female
  out <- bootmakr(f, data = darfur_data, treat = "directlyharmed",
                  benchmark_covariates = "female", reps = 30, seed = 1, progress = FALSE)
  kc  <- out$breakdown$kd_conventional
  ob  <- sensemakr::ovb_bounds(lm(f, data = darfur_data), treatment = "directlyharmed",
                               benchmark_covariates = "female",
                               kd = kc * c(0.5, 0.99, 1 - 1e-7, 1 + 1e-7))
  outside <- ob$adjusted_lower_CI * ob$adjusted_upper_CI > 0
  expect_identical(outside, c(TRUE, TRUE, TRUE, FALSE))
})

test_that("it works with weights, a grouped benchmark and clusters", {
  d <- firms
  d$w  <- 0.5 + (d$year - min(d$year)) / diff(range(d$year))
  args <- list(y ~ x + c, data = d, treat = "x", gbenchmark_covariates = "c",
               weights = "w", cluster = "firm", reps = 60, seed = 5)
  out <- do.call(bootmakr, c(args, list(progress = FALSE)))
  bd  <- out$breakdown$kd

  expect_identical(out$breakdown$benchmark, "c")
  expect_lt(do.call(p_at, c(list(bd * (1 - 1e-8)), args)), 0.05)
  expect_gte(do.call(p_at, c(list(bd * (1 + 1e-8)), args)), 0.05)
})

test_that("with two-way clustering it is where the z-test stops rejecting", {
  args <- list(y ~ x + c, data = firms, treat = "x", benchmark_covariates = "c",
               cluster = c("firm", "year"), reps = 40, seed = 6)
  out <- do.call(bootmakr, c(args, list(progress = FALSE)))
  bd  <- out$breakdown

  expect_identical(bd$method, "cgm_twoway")
  expect_identical(names(out$boot_bench), c("G", "H", "GH"))
  expect_true(bd$kd > 0 && bd$kd < bd$kd_limit)
  expect_lt(do.call(p_at, c(list(bd$kd * (1 - 1e-4)), args)), 0.05)
  expect_gt(do.call(p_at, c(list(bd$kd * (1 + 1e-4)), args)), 0.05)
})

test_that("replicates beyond their own limit are left out, as in the table", {
  d    <- coll_data()
  args <- list(y ~ x + c, data = d, treat = "x", benchmark_covariates = "c",
               reps = 100, seed = 11)
  out <- suppressWarnings(do.call(bootmakr, c(args, list(kd = 0.1, progress = FALSE))))
  bd  <- out$breakdown
  left_out <- bd$n_reps - bd$n_used

  expect_gt(left_out, 0)
  expect_match(bd$note, "left out")
  lo <- suppressWarnings(do.call(bootmakr, c(args, list(kd = bd$kd * (1 - 1e-8), progress = FALSE))))
  hi <- suppressWarnings(do.call(bootmakr, c(args, list(kd = bd$kd * (1 + 1e-8), progress = FALSE))))
  expect_lt(lo$results$pvalue, 0.05)
  expect_gte(hi$results$pvalue, 0.05)
  expect_equal(hi$N_fail, left_out)
})

test_that("replicates with a coefficient of the opposite sign are handled", {
  set.seed(11)
  d <- data.frame(x = rnorm(150), c = rnorm(150))
  d$y <- 0.10 * d$x + 0.5 * d$c + rnorm(150)
  args <- list(y ~ x + c, data = d, treat = "x", benchmark_covariates = "c",
               alpha = 0.30, reps = 60, seed = 11)
  out <- do.call(bootmakr, c(args, list(progress = FALSE)))
  bd  <- out$breakdown$kd

  expect_gt(sum(out$boot_base[, "estimate"] < 0), 0)
  expect_gt(bd, 0)
  expect_lt(do.call(p_at, c(list(bd * (1 - 1e-8)), args)), 0.30)
  expect_gte(do.call(p_at, c(list(bd * (1 + 1e-8)), args)), 0.30)
})

test_that("cases without a breakdown kd say why", {
  skip_helper()
  f <- peacefactor ~ directlyharmed + age + female
  manual <- bootmakr(f, data = darfur_data, treat = "directlyharmed",
                     r2dz.x = 0.02, r2yz.dx = 0.05, reps = 20, seed = 7, progress = FALSE)
  expect_true(is.na(manual$breakdown$kd))
  expect_true(is.na(manual$breakdown$benchmark))
  expect_true(all(is.na(manual$boot_bench)))
  expect_false(any(grepl("Breakdown kd", capture.output(print(manual)), fixed = TRUE)))
  # the robustness value does not need a benchmark
  expect_false(is.na(manual$robustness$rv_qa_boot))

  up <- bootmakr(f, data = darfur_data, treat = "directlyharmed",
                 benchmark_covariates = "female", reduce = FALSE,
                 reps = 20, seed = 7, progress = FALSE)
  expect_true(is.na(up$breakdown$kd))
  expect_match(up$breakdown$note, "reduce = TRUE")

  other_ky <- bootmakr(f, data = darfur_data, treat = "directlyharmed",
                       benchmark_covariates = "female", kd = 1, ky = 2,
                       reps = 20, seed = 7, progress = FALSE)
  expect_false(is.na(other_ky$breakdown$kd))
  expect_match(other_ky$breakdown$note, "ky = kd")

  # second benchmark selected with bounds_row
  second <- bootmakr(f, data = darfur_data, treat = "directlyharmed",
                     benchmark_covariates = c("age", "female"), bounds_row = 2,
                     reps = 20, seed = 7, progress = FALSE)
  first  <- bootmakr(f, data = darfur_data, treat = "directlyharmed",
                     benchmark_covariates = "female", reps = 20, seed = 7, progress = FALSE)
  expect_identical(second$breakdown$benchmark, "female")
  expect_equal(second$breakdown$kd, first$breakdown$kd, tolerance = 1e-12)

  # an estimate that is not significant to begin with breaks down at zero
  set.seed(5)
  d <- data.frame(x = rnorm(150), c = rnorm(150))
  d$y <- 0.5 * d$c + rnorm(150)
  ns <- bootmakr(y ~ x + c, data = d, treat = "x", benchmark_covariates = "c",
                 reps = 60, seed = 8, progress = FALSE)
  expect_identical(ns$breakdown$kd, 0)
})

test_that("the search reports when nothing crosses or the replicates run out", {
  base <- cbind(estimate = rep(1, 50), scale = rep(1, 50))
  # benchmark unrelated to the outcome: no bias at any kd
  none <- bootmakr:::.breakdown_percentile(
    base, cbind(r2dxj.x = rep(0.01, 50), r2yxj.dx = rep(0, 50)), 0.05, k_limit = 50)
  expect_true(is.na(none$kd))
  expect_match(none$note, "stays significant")
  # replicates reach their limits one after the other without crossing
  out <- bootmakr:::.breakdown_percentile(
    base, cbind(r2dxj.x = seq(0.3, 0.6, length.out = 50), r2yxj.dx = rep(0, 50)), 0.05, k_limit = 50)
  expect_true(is.na(out$kd))
  expect_match(out$note, "Too few")
})

test_that("print shows the breakdown kd", {
  skip_helper()
  out <- bootmakr(peacefactor ~ directlyharmed + age + female, data = darfur_data,
                  treat = "directlyharmed", benchmark_covariates = "female",
                  reps = 30, seed = 9, progress = FALSE)
  txt <- capture.output(print(out))
  expect_true(any(grepl("Breakdown kd", txt, fixed = TRUE)))
  expect_true(any(grepl("largest kd this benchmark allows", txt, fixed = TRUE)))
})
