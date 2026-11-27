# Semi-weights for the benchmark (Wainstein & Hazlett, 2026): the benchmark's
# partial R2 with the treatment is measured under the weights and under the
# semi-weights, and the bound uses both.

data("firms", package = "bootmakr")
firms$w  <- exp(0.4 * as.numeric(scale(firms$c)))     # weights that depend on the benchmark
firms$s  <- 1                                         # ... and semi-weights formed without it
firms$s2 <- 0.5 + (firms$year - min(firms$year)) / diff(range(firms$year))
cols <- c("estimate", "se", "ci_lower", "ci_upper", "pvalue")
run  <- function(...) suppressWarnings(bootmakr(y ~ x + c, firms, treat = "x", ..., progress = FALSE))

# the bound of Wainstein and Hazlett (2026, equations 25 and 26), by hand
wh_bound <- function(r2w, r2s, r2y, kd, ky = kd) {
  r2dz <- kd * r2s / (1 - r2w)
  r2zx <- kd * r2s / (1 - kd * r2s) * r2w / (1 - r2w)
  list(r2dz.x = r2dz,
       r2yz.dx = ((sqrt(ky) + sqrt(r2zx)) / sqrt(1 - r2zx))^2 * r2y / (1 - r2y))
}

test_that("the bound with semi-weights is equations (25) and (26), computed by hand", {
  out <- run(benchmark_covariates = "c", kd = c(1, 2), ky = c(2, 1), weights = "w",
             semiweights = "s", reps = 30, seed = 1)
  fit_y <- lm(y ~ x + c, firms, weights = w)
  r2w <- unname(sensemakr::partial_r2(lm(x ~ c, firms, weights = w), covariates = "c"))
  r2s <- unname(sensemakr::partial_r2(lm(x ~ c, firms, weights = s), covariates = "c"))
  r2y <- unname(sensemakr::partial_r2(fit_y, covariates = "c"))
  ob  <- out$benchmark_strength$observed
  expect_equal(c(ob$r2dxj.x, ob$r2dxj.x_semi, ob$r2yxj.dx), c(r2w, r2s, r2y), tolerance = 1e-10)
  expect_equal(ob$r_dxj.x_semi, sqrt(r2s), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(r2w, r2s)))

  bnd <- wh_bound(r2w, r2s, r2y, kd = c(1, 2), ky = c(2, 1))
  expect_equal(out$results$r2dz.x,  bnd$r2dz.x,  tolerance = 1e-10)
  expect_equal(out$results$r2yz.dx, bnd$r2yz.dx, tolerance = 1e-10)
  hand <- sensemakr::adjusted_estimate(fit_y, treatment = "x", r2dz.x = bnd$r2dz.x,
                                       r2yz.dx = bnd$r2yz.dx)
  expect_equal(out$results$estimate, unname(hand), tolerance = 1e-10)
  expect_true(out$semiweights)
  expect_output(print(out), "semi-weights")
  expect_output(print(out), "with treatment \\| X, semi-weights")
})

test_that("semi-weights equal to the weights give the results without semi-weights", {
  for (cl in list(NULL, "firm", c("firm", "year"))) {
    a <- run(benchmark_covariates = "c", kd = c(1, 1.5), cluster = cl, weights = "w", reps = 40, seed = 2)
    b <- run(benchmark_covariates = "c", kd = c(1, 1.5), cluster = cl, weights = "w", semiweights = "w",
             reps = 40, seed = 2)
    expect_equal(b$results, a$results, tolerance = 1e-12)
    expect_equal(b$boot_samples, a$boot_samples, tolerance = 1e-12)
    expect_equal(b$breakdown[c("kd", "kd_conventional", "kd_limit")],
                 a$breakdown[c("kd", "kd_conventional", "kd_limit")], tolerance = 1e-10)
    expect_false(a$semiweights)
  }
  # the stored benchmark R2 gain a third column, equal to the first here
  expect_identical(colnames(b$boot_bench$G), c("r2dxj.x", "r2yxj.dx", "r2dxj.x_semi"))
  expect_equal(b$boot_bench$G[, 3], b$boot_bench$G[, 1], tolerance = 1e-12)
  expect_identical(colnames(a$boot_bench$G), c("r2dxj.x", "r2yxj.dx"))
})

test_that("semi-weights move the bound only: coefficient, scale and robustness value stay", {
  a <- run(benchmark_covariates = "c", cluster = "firm", weights = "w", reps = 40, seed = 3)
  b <- run(benchmark_covariates = "c", cluster = "firm", weights = "w", semiweights = "s", reps = 40, seed = 3)
  expect_identical(b$boot_base, a$boot_base)
  expect_identical(b$robustness, a$robustness)
  expect_false(isTRUE(all.equal(b$results$estimate, a$results$estimate)))
  expect_equal(dim(b$boot_bench), c(40L, 3L))
  # every replicate: adjusted estimate from its own three partial R2
  bnd <- wh_bound(b$boot_bench[, 1], b$boot_bench[, 3], b$boot_bench[, 2], kd = 1)
  adj <- sign(b$boot_base[, 1]) * (abs(b$boot_base[, 1]) -
    b$boot_base[, 2] * sqrt(bnd$r2yz.dx * bnd$r2dz.x / (1 - bnd$r2dz.x)))
  expect_equal(unname(b$boot_samples[, 1]), unname(adj), tolerance = 1e-10)
})

test_that("semi-weight functions are called on every resample and equal stored columns", {
  seen <- integer(0)
  spy  <- function(d) { seen <<- c(seen, nrow(d)); d$s2 }
  a <- run(benchmark_covariates = "c", cluster = "firm", weights = "w", semiweights = "s2", reps = 30, seed = 4)
  b <- run(benchmark_covariates = "c", cluster = "firm", weights = function(d) d$w, semiweights = spy,
           reps = 30, seed = 4)
  expect_length(seen, 31L)                      # original data and every resample
  expect_equal(b$results, a$results)
  expect_equal(b$boot_samples, a$boot_samples)
  expect_equal(b$boot_bench, a$boot_bench)
  expect_equal(b$breakdown, a$breakdown)

  # semi-weights that depend on the sample they are computed on
  sfun <- function(d) 1 / (1 + (d$x - mean(d$x))^2)
  f <- run(benchmark_covariates = "c", cluster = "firm", weights = function(d) d$w, semiweights = sfun,
           reps = 30, seed = 4)
  g <- run(benchmark_covariates = "c", cluster = "firm", weights = "w", semiweights = sfun(firms),
           reps = 30, seed = 4)
  expect_equal(f$results$estimate, g$results$estimate)
  expect_identical(f$boot_base, g$boot_base)
  expect_false(isTRUE(all.equal(f$boot_bench[, 3], g$boot_bench[, 3])))
})

test_that("fix_benchmark with semi-weights equals a manual-R2 run at the original bound", {
  d0 <- run(benchmark_covariates = "c", cluster = "firm", weights = "w", semiweights = "s", reps = 50, seed = 5)
  fx <- run(benchmark_covariates = "c", cluster = "firm", weights = "w", semiweights = "s", reps = 50, seed = 5,
            fix_benchmark = TRUE)
  mn <- run(r2dz.x = d0$results$r2dz.x, r2yz.dx = d0$results$r2yz.dx, cluster = "firm", weights = "w",
            reps = 50, seed = 5)
  expect_equal(fx$results[cols], mn$results[cols], tolerance = 1e-10)
  expect_equal(unname(fx$boot_samples), unname(mn$boot_samples), tolerance = 1e-10)
  expect_true(all(fx$boot_bench[, 3] == fx$boot_bench[1, 3]))
  expect_equal(fx$boot_bench[1, 3], d0$benchmark_strength$observed$r2dxj.x_semi, ignore_attr = TRUE)
  # a semi-weight function is not called again when the benchmark is fixed
  n <- 0
  fxf <- run(benchmark_covariates = "c", cluster = "firm", weights = function(d) d$w,
             semiweights = function(d) { n <<- n + 1; d$s }, reps = 50, seed = 5, fix_benchmark = TRUE)
  expect_equal(n, 1)
  expect_equal(fxf$results[cols], fx$results[cols], tolerance = 1e-10)
})

test_that("the breakdown kd holds with semi-weights, and the limit is (1 - R2w) / R2s", {
  for (cl in list("firm", c("firm", "year"))) {
    args <- list(benchmark_covariates = "c", cluster = cl, weights = "w", semiweights = "s", reps = 100, seed = 6)
    out <- do.call(run, args)
    k   <- out$breakdown$kd
    expect_true(is.finite(k) && k > 0)
    expect_lt(do.call(run, c(args, list(kd = k * (1 - 1e-6))))$results$pvalue, 0.05)
    expect_gte(do.call(run, c(args, list(kd = k * (1 + 1e-6))))$results$pvalue, 0.05)
  }
  ob  <- out$benchmark_strength$observed
  lim <- (1 - ob$r2dxj.x) / ob$r2dxj.x_semi
  expect_equal(out$breakdown$kd_limit, lim, tolerance = 1e-12)
  expect_error(bootmakr(y ~ x + c, firms, treat = "x", benchmark_covariates = "c", weights = "w",
                        semiweights = "s", kd = lim * 1.001, reps = 20, progress = FALSE),
               "Impossible kd value")
  expect_s3_class(run(benchmark_covariates = "c", weights = "w", semiweights = "s", kd = lim * 0.999,
                      reps = 20, seed = 1), "bootmakr")
})

test_that("a covariate gives the same answer as a single and as a grouped benchmark", {
  for (cl in list(NULL, "firm")) {
    single <- run(benchmark_covariates = "c", kd = c(1, 2), cluster = cl, weights = "w", semiweights = "s2",
                  reps = 30, seed = 7)
    group  <- run(gbenchmark_covariates = "c", kd = c(1, 2), cluster = cl, weights = "w", semiweights = "s2",
                  reps = 30, seed = 7)
    expect_equal(single$results$estimate, group$results$estimate, tolerance = 1e-10)
    expect_equal(single$results$r2dz.x,   group$results$r2dz.x,   tolerance = 1e-10)
    expect_equal(single$results$r2yz.dx,  group$results$r2yz.dx,  tolerance = 1e-10)
    expect_equal(unname(single$boot_samples), unname(group$boot_samples), tolerance = 1e-10)
    expect_equal(unname(single$boot_bench), unname(group$boot_bench), tolerance = 1e-10)
  }
})

test_that("several benchmarks take a list with one set of semi-weights each", {
  skip_helper()
  fml <- peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female
  d    <- darfur_data
  d$w  <- ifelse(d$female == 1, 2, 0.5) * ifelse(d$pastvoted == 1, 1.5, 1)
  d$sf <- ifelse(d$pastvoted == 1, 1.5, 1)      # formed without female
  d$sp <- ifelse(d$female == 1, 2, 0.5)         # formed without pastvoted
  one <- function(...) suppressWarnings(bootmakr(fml, d, treat = "directlyharmed", weights = "w",
                                                 reps = 30, seed = 8, progress = FALSE, ...))
  both <- one(benchmark_covariates = c("female", "pastvoted"), semiweights = list("sf", d$sp), bounds_row = 2)
  nm   <- one(benchmark_covariates = c("female", "pastvoted"), semiweights = c("sf", "sp"), bounds_row = 2)
  pv   <- one(benchmark_covariates = "pastvoted", semiweights = "sp")
  fe   <- one(benchmark_covariates = "female", semiweights = "sf")
  expect_equal(both$results, nm$results)
  expect_equal(both$results$estimate, pv$results$estimate, tolerance = 1e-10)
  expect_equal(both$results$r2dz.x,   pv$results$r2dz.x,   tolerance = 1e-10)
  expect_equal(both$benchmark_strength$observed$r2dxj.x_semi,
               c(fe$benchmark_strength$observed$r2dxj.x_semi, pv$benchmark_strength$observed$r2dxj.x_semi),
               tolerance = 1e-10)
  expect_equal(both$breakdown$kd, pv$breakdown$kd, tolerance = 1e-8)
})

test_that("rows dropped for missing values keep weights and semi-weights aligned", {
  f2 <- firms; f2$y[c(3, 50, 200)] <- NA; f2$c[c(10, 400)] <- NA
  out <- suppressWarnings(bootmakr(y ~ x + c, f2, treat = "x", benchmark_covariates = "c", weights = "w",
                                   semiweights = "s2", reps = 30, seed = 9, progress = FALSE))
  cc  <- complete.cases(f2[, c("y", "x", "c")])
  r2s <- sensemakr::partial_r2(lm(x ~ c, f2[cc, ], weights = s2), covariates = "c")
  expect_equal(out$benchmark_strength$observed$r2dxj.x_semi, unname(r2s), tolerance = 1e-10)
  expect_equal(out$N, sum(cc))
  expect_equal(out$N_fail, 0L)
})

test_that("misuse of semiweights stops; a failing function is a failed replication", {
  expect_error(run(benchmark_covariates = "c", semiweights = "s", reps = 20), "needs `weights`")
  expect_error(run(r2dz.x = 0.01, r2yz.dx = 0.1, weights = "w", semiweights = "s", reps = 20),
               "needs a benchmark")
  expect_error(run(benchmark_covariates = "c", weights = "w", semiweights = list("s", "s2"), reps = 20),
               "one entry per benchmark covariate")
  expect_error(run(benchmark_covariates = "c", weights = "w", semiweights = function(d) d$s, reps = 20),
               "must be one too")
  expect_error(run(benchmark_covariates = "c", weights = function(d) d$w, semiweights = "s", reps = 20),
               "must be a function too")
  expect_error(run(benchmark_covariates = "c", weights = "w", semiweights = -firms$s2, reps = 20),
               "non-negative")
  expect_error(run(benchmark_covariates = "c", weights = "w", semiweights = replace(firms$s2, 4, NA), reps = 20),
               "non-negative")
  expect_error(run(benchmark_covariates = "c", weights = function(d) d$w, semiweights = function(d) d$s[-1],
                   reps = 20), "`semiweights` function must return")

  n <- 0
  flaky <- function(d) { n <<- n + 1; if (n %% 5 == 0) stop("no semi-weights"); d$s }
  out <- run(benchmark_covariates = "c", weights = function(d) d$w, semiweights = flaky, reps = 40, seed = 10)
  expect_equal(out$N_fail, 8L)
  expect_equal(sum(!is.finite(out$boot_samples[, 1])), 8L)
  expect_true(all(is.finite(out$boot_base)))       # coefficient and scale are kept
})
