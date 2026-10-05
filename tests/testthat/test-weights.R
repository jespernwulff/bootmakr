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

  # same adjusted estimate as sensemakr on the weighted regression
  fit_w <- lm(fml, data = d, weights = w)
  expect_equal(
    out_v$results$estimate,
    sensemakr::ovb_bounds(fit_w, treatment = "directlyharmed",
                          benchmark_covariates = "female",
                          kd = 1)$adjusted_estimate,
    tolerance = 1e-10
  )
  expect_false(isTRUE(all.equal(out_v$results$estimate, out_u$results$estimate)))

  # no replication is lost, and both ways of passing weights agree
  expect_equal(out_v$N_fail, 0L)
  expect_equal(out_v$N_successful, 40L)
  expect_identical(out_v$results, out_c$results)
  expect_identical(out_v$boot_samples, out_c$boot_samples)
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
