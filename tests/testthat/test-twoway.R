test_that("two-way cluster bootstrap returns CGM structure", {
  skip_helper()

  out <- suppressWarnings(bootmakr(
    darfur_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    benchmark_covariates = "female",
    cluster = c("village", "female"),
    reps    = 100,
    seed    = 3,
    progress = FALSE
  ))

  expect_identical(out$method, "cgm_twoway")
  expect_true(is.list(out$boot_samples))
  expect_named(out$boot_samples, c("G", "H", "GH"))
  expect_true(all(vapply(out$boot_samples, is.matrix, logical(1))))
  expect_equal(nrow(out$boot_samples$G), 100L)
  expect_equal(nrow(out$boot_samples$H), 100L)
  expect_equal(nrow(out$boot_samples$GH), 100L)
  expect_equal(unname(out$N_clust["G"]),  486L)
  expect_equal(unname(out$N_clust["H"]),  2L)
  expect_equal(unname(out$N_clust["GH"]), 576L)
  expect_identical(out$cluster_names, c("village", "female"))
  # var_G/var_H/var_GH columns exist
  expect_true(all(c("var_G", "var_H", "var_GH", "var_neg_fix") %in%
                    names(out$results)))
  # Point estimate is unchanged from one-way (it's the original fit)
  expect_equal(out$results$estimate, 0.07522, tolerance = 1e-4)
  # SE is finite and positive
  expect_true(is.finite(out$results$se))
  expect_true(out$results$se > 0)
})

test_that("two-way SE matches V_G + V_H - V_GH (or max() fallback)", {
  skip_helper()
  out <- suppressWarnings(bootmakr(
    darfur_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    benchmark_covariates = "female",
    cluster = c("village", "female"),
    reps    = 100,
    seed    = 4,
    progress = FALSE
  ))
  r <- out$results
  v_raw <- r$var_G + r$var_H - r$var_GH
  expected_var <- if (isTRUE(r$var_neg_fix)) max(r$var_G, r$var_H) else v_raw
  expect_equal(r$se^2, expected_var, tolerance = 1e-10)
})

test_that("two-way emits soft warning when a dimension has < 10 clusters", {
  skip_helper()
  # Collect all warnings (the negative-variance fallback may also fire at
  # small reps, so a plain expect_warning would only catch the first).
  warnings_seen <- character()
  res <- withCallingHandlers(
    bootmakr(
      darfur_formula,
      data    = darfur_data,
      treat   = "directlyharmed",
      benchmark_covariates = "female",
      cluster = c("village", "female"),
      reps    = 50,
      seed    = 5,
      progress = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("only 2 unique", warnings_seen)),
              info = paste("Got warnings:",
                           paste(warnings_seen, collapse = " | ")))
})

test_that("two-way kd sweep produces one row per kd", {
  skip_helper()
  out <- suppressWarnings(bootmakr(
    darfur_formula,
    data    = darfur_data,
    treat   = "directlyharmed",
    benchmark_covariates = "female",
    kd      = c(0.5, 1, 2),
    cluster = c("village", "female"),
    reps    = 60,
    seed    = 6,
    progress = FALSE
  ))
  expect_equal(nrow(out$results), 3L)
  expect_equal(out$results$kd, c(0.5, 1, 2))
})
