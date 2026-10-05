test_that("one-way cluster bootstrap returns expected object structure", {
  skip_helper()

  out <- bootmakr(
    darfur_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    benchmark_covariates = "female",
    cluster = "village",
    reps    = 100,
    seed    = 1,
    progress = FALSE
  )

  expect_s3_class(out, "bootmakr")
  expect_identical(out$method, "percentile")
  expect_true(is.matrix(out$boot_samples))
  expect_equal(nrow(out$boot_samples), 100L)
  expect_equal(ncol(out$boot_samples), 1L)
  expect_equal(out$N_clust, 486L)
  expect_identical(out$cluster_names, "village")
  expect_equal(nrow(out$results), 1L)
  expect_true(all(c("kd", "estimate", "se", "ci_lower",
                    "ci_upper", "pvalue") %in% names(out$results)))
  expect_equal(out$results$estimate, 0.07522, tolerance = 1e-4)
  expect_true(out$results$se > 0)
})

test_that("kd sweep yields one row per kd", {
  skip_helper()
  out <- bootmakr(
    darfur_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    benchmark_covariates = "female",
    kd      = c(0.5, 1, 1.5),
    cluster = "village",
    reps    = 80,
    seed    = 2,
    progress = FALSE
  )
  expect_equal(nrow(out$results), 3L)
  expect_equal(out$results$kd, c(0.5, 1, 1.5))
})

test_that("replication counts of 10,000 and more print as integers", {
  skip_helper()
  out <- bootmakr(
    peacefactor ~ directlyharmed + female,
    data    = darfur_data,
    treat   = "directlyharmed",
    benchmark_covariates = "female",
    reps    = 20,
    seed    = 7,
    progress = FALSE
  )
  out$N_reps <- 10000
  txt <- capture.output(print(out))
  expect_true(any(grepl("(10,000 reps", txt, fixed = TRUE)))
  expect_false(any(grepl("1e+04", txt, fixed = TRUE)))
})
