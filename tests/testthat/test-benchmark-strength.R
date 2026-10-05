# The benchmark-strength block is descriptive: it must agree with sensemakr's
# own benchmarking quantities and must not depend on the bootstrap.

light_formula <- peacefactor ~ directlyharmed + age + farmer_dar +
  herder_dar + pastvoted + hhsize_darfur + female

test_that("single benchmark: observed and implied strength agree with sensemakr", {
  skip_helper()
  out <- bootmakr(
    light_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    benchmark_covariates = "female",
    kd      = c(0.5, 1, 2),
    cluster = "village",
    reps    = 60,
    seed    = 11,
    progress = FALSE
  )
  bs <- out$benchmark_strength
  expect_identical(bs$type, "single")

  fit   <- lm(light_formula, data = darfur_data)
  fit_d <- lm(update(light_formula, directlyharmed ~ . - directlyharmed),
              data = darfur_data)

  # observed strength of the benchmark
  expect_equal(bs$observed$r2yxj.dx,
               unname(sensemakr::partial_r2(fit, covariates = "female")),
               tolerance = 1e-10)
  expect_equal(bs$observed$r2dxj.x,
               unname(sensemakr::partial_r2(fit_d, covariates = "female")),
               tolerance = 1e-10)
  expect_equal(bs$observed$r_dxj.x,  sqrt(bs$observed$r2dxj.x))
  expect_equal(bs$observed$r_yxj.dx, sqrt(bs$observed$r2yxj.dx))

  # implied strength of the omitted variable = sensemakr's bounds
  bnds <- sensemakr::ovb_bounds(fit, treatment = "directlyharmed",
                                benchmark_covariates = "female",
                                kd = c(0.5, 1, 2), ky = c(0.5, 1, 2))
  expect_equal(bs$implied$kd, c(0.5, 1, 2))
  expect_equal(bs$implied$r2dz.x,  bnds$r2dz.x,  tolerance = 1e-10)
  expect_equal(bs$implied$r2yz.dx, bnds$r2yz.dx, tolerance = 1e-10)
  expect_equal(bs$implied$r_dz.x,  sqrt(bs$implied$r2dz.x))
  expect_equal(bs$implied$r_yz.dx, sqrt(bs$implied$r2yz.dx))

  # r2dz.x = kd * r2dxj.x / (1 - r2dxj.x) (Cinelli & Hazlett, 2020)
  expect_equal(bs$implied$r2dz.x,
               c(0.5, 1, 2) * bs$observed$r2dxj.x / (1 - bs$observed$r2dxj.x),
               tolerance = 1e-10)

  # the results table carries the implied strength each row is based on,
  # and the adjusted estimates are sensemakr's
  expect_equal(out$results$r2dz.x,  bs$implied$r2dz.x)
  expect_equal(out$results$r2yz.dx, bs$implied$r2yz.dx)
  expect_equal(out$results$estimate, bnds$adjusted_estimate, tolerance = 1e-10)
})

test_that("benchmark strength does not depend on seed, reps or resampling scheme", {
  skip_helper()
  args <- list(light_formula, data = darfur_data, treat = "directlyharmed",
               benchmark_covariates = "female", progress = FALSE)
  a <- do.call(bootmakr, c(args, list(reps = 40, seed = 1)))
  b <- do.call(bootmakr, c(args, list(reps = 60, seed = 2, cluster = "village")))
  expect_identical(a$benchmark_strength, b$benchmark_strength)
})

test_that("grouped benchmark: strength uses the group partial R2", {
  skip_helper()
  grp <- c("female", "pastvoted")
  out <- bootmakr(
    light_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    gbenchmark_covariates = grp,
    kd      = c(0.5, 1),
    reps    = 40,
    seed    = 12,
    progress = FALSE
  )
  bs <- out$benchmark_strength
  expect_identical(bs$type, "group")
  expect_equal(nrow(bs$observed), 1L)
  expect_identical(bs$observed$benchmark, "female + pastvoted")

  fit   <- lm(light_formula, data = darfur_data)
  fit_d <- lm(update(light_formula, directlyharmed ~ . - directlyharmed),
              data = darfur_data)
  expect_equal(bs$observed$r2yxj.dx,
               unname(sensemakr::group_partial_r2(fit, covariates = grp)),
               tolerance = 1e-10)
  expect_equal(bs$observed$r2dxj.x,
               unname(sensemakr::group_partial_r2(fit_d, covariates = grp)),
               tolerance = 1e-10)
  expect_equal(nrow(bs$implied), 2L)
  expect_equal(out$results$r2dz.x, bs$implied$r2dz.x)
})

test_that("manual R2: implied strength is the supplied pair, no observed block", {
  skip_helper()
  out <- bootmakr(
    light_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    r2dz.x  = 0.01,
    r2yz.dx = 0.02,
    reps    = 40,
    seed    = 13,
    progress = FALSE
  )
  bs <- out$benchmark_strength
  expect_identical(bs$type, "manual")
  expect_null(bs$observed)
  expect_equal(bs$implied$r2dz.x,  0.01)
  expect_equal(bs$implied$r2yz.dx, 0.02)
  expect_equal(out$results$r2dz.x, 0.01)
})

test_that("print shows the benchmark-strength block", {
  skip_helper()
  out <- bootmakr(
    light_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    benchmark_covariates = "female",
    reps    = 40,
    seed    = 14,
    progress = FALSE
  )
  txt <- capture.output(print(out))
  expect_true(any(grepl("^Benchmark strength", txt)))
  expect_true(any(grepl("female with treatment \\| X", txt)))
  expect_true(any(grepl("Implied strength of the omitted variable", txt)))
})
