test_that("two-way clustering with strata errors", {
  skip_helper()
  expect_error(
    bootmakr(
      peacefactor ~ directlyharmed + female,
      data    = darfur_data,
      treat   = "directlyharmed",
      benchmark_covariates = "female",
      cluster = c("village", "female"),
      strata  = "female",
      reps    = 10,
      progress = FALSE
    ),
    regexp = "strata"
  )
})

test_that("cluster length > 2 errors", {
  skip_helper()
  expect_error(
    bootmakr(
      peacefactor ~ directlyharmed + female,
      data    = darfur_data,
      treat   = "directlyharmed",
      benchmark_covariates = "female",
      cluster = c("village", "female", "age"),
      reps    = 10,
      progress = FALSE
    ),
    regexp = "two"
  )
})

test_that("missing benchmark spec errors", {
  skip_helper()
  expect_error(
    bootmakr(
      peacefactor ~ directlyharmed + female,
      data    = darfur_data,
      treat   = "directlyharmed",
      reps    = 10,
      progress = FALSE
    ),
    regexp = "benchmark_covariates|r2dz|r2yz|gbenchmark"
  )
})

test_that("unknown cluster column errors", {
  skip_helper()
  expect_error(
    bootmakr(
      peacefactor ~ directlyharmed + female,
      data    = darfur_data,
      treat   = "directlyharmed",
      benchmark_covariates = "female",
      cluster = "no_such_column",
      reps    = 10,
      progress = FALSE
    ),
    regexp = "not found"
  )
})
