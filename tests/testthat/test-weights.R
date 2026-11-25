# `weights` must reach lm() whether it is supplied as a vector or as a column
# name, and whatever environment the formula was created in.

test_that("weights are used in the original fit and in every replication", {
  skip_helper()
  fml <- peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female
  d   <- darfur_data
  d$w <- seq(0.5, 1.5, length.out = nrow(d))

  run <- function(...) bootmakr(fml, data = d, treat = "directlyharmed",
                                benchmark_covariates = "female",
                                reps = 40, seed = 21, progress = FALSE, ...)
  out_u <- run()
  out_v <- run(weights = d$w)     # vector
  out_c <- run(weights = "w")     # column name

  # The adjusted estimate is sensemakr's bound formula applied to partial R2s
  # that are weighted on both sides: the outcome model and the treatment
  # model each carry the weights.
  fit_w  <- lm(fml, data = d, weights = w)
  fit_dw <- lm(update(fml, directlyharmed ~ . - directlyharmed), data = d, weights = w)
  bnd <- sensemakr::ovb_partial_r2_bound(
    r2dxj.x  = sensemakr::partial_r2(fit_dw, covariates = "female"),
    r2yxj.dx = sensemakr::partial_r2(fit_w,  covariates = "female"),
    kd = 1
  )
  hand <- sensemakr::adjusted_estimate(fit_w, treatment = "directlyharmed",
                                       r2dz.x = bnd$r2dz.x, r2yz.dx = bnd$r2yz.dx)
  expect_equal(out_v$results$estimate, unname(hand), tolerance = 1e-10)
  expect_equal(out_v$results$r2dz.x, unname(bnd$r2dz.x), tolerance = 1e-10)
  expect_equal(out_v$benchmark_strength$observed$r2dxj.x,
               unname(sensemakr::partial_r2(fit_dw, covariates = "female")),
               tolerance = 1e-10)

  expect_false(isTRUE(all.equal(out_v$results$estimate, out_u$results$estimate)))

  # no replication is lost, and both ways of passing weights agree
  expect_equal(out_v$N_fail, 0L)
  expect_equal(out_v$N_successful, 40L)
  expect_identical(out_v$results, out_c$results)
  expect_identical(out_v$boot_samples, out_c$boot_samples)
})

test_that("a covariate gives the same answer as a single and as a grouped benchmark", {
  skip_helper()
  fml <- peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female
  d    <- darfur_data
  d$w  <- seq(0.5, 1.5, length.out = nrow(d))
  d$wf <- ifelse(d$female == 1, 2, 0.5)   # weights that depend on the benchmark

  for (wn in list(NULL, "w", "wf")) {
    run <- function(...) bootmakr(fml, data = d, treat = "directlyharmed",
                                  kd = c(1, 2), weights = wn,
                                  reps = 30, seed = 23, progress = FALSE, ...)
    single <- run(benchmark_covariates = "female")
    group  <- run(gbenchmark_covariates = "female")
    expect_equal(single$results$estimate, group$results$estimate, tolerance = 1e-10)
    expect_equal(single$results$r2dz.x,   group$results$r2dz.x,   tolerance = 1e-10)
    expect_equal(single$results$r2yz.dx,  group$results$r2yz.dx,  tolerance = 1e-10)
    expect_equal(single$results$impact,   group$results$impact,   tolerance = 1e-10)
    expect_equal(unname(single$boot_samples), unname(group$boot_samples),
                 tolerance = 1e-10)
  }
})

test_that("without weights the adjusted estimate is sensemakr's", {
  skip_helper()
  fml <- peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female
  out <- bootmakr(fml, data = darfur_data, treat = "directlyharmed",
                  benchmark_covariates = "female", kd = c(1, 2, 3),
                  reps = 20, seed = 24, progress = FALSE)
  sm <- sensemakr::ovb_bounds(lm(fml, data = darfur_data),
                              treatment = "directlyharmed",
                              benchmark_covariates = "female", kd = c(1, 2, 3))
  expect_equal(out$results$estimate, sm$adjusted_estimate, tolerance = 1e-10)
  expect_equal(out$results$r2dz.x,   sm$r2dz.x,            tolerance = 1e-10)
  expect_equal(out$results$r2yz.dx,  sm$r2yz.dx,           tolerance = 1e-10)
})

test_that("weights work with a grouped benchmark and a cluster bootstrap", {
  skip_helper()
  fml <- peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
    pastvoted + hhsize_darfur + female
  w <- seq(0.5, 1.5, length.out = nrow(darfur_data))
  out <- bootmakr(fml, data = darfur_data, treat = "directlyharmed",
                  gbenchmark_covariates = c("female", "pastvoted"),
                  weights = w, cluster = "village",
                  reps = 40, seed = 22, progress = FALSE)
  expect_equal(out$N_fail, 0L)
  expect_true(is.finite(out$results$se) && out$results$se > 0)
})
